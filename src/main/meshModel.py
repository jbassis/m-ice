from   fenics import *
import numpy  as np
from mshr import *
import os
import math
import pygmsh
#Meshio for saving mesh to file for use in fenics
import meshio

from importlib import reload

import mesh_functions
reload(mesh_functions)
from mesh_functions import *

import os
import tempfile


def meshGmsh(geometryArray, meshSize):

    #Define geometry class which constains all functions to add different features (points, lines, surfaces, etc)
    geom = pygmsh.built_in.geometry.Geometry()

    #Iterate through each coordinate in the geometryArray, making a point with each (assumed constant mesh size) and adding it to the points array.
    points = []
    for j in range(np.shape(geometryArray)[1]):
        points.append(geom.add_point([geometryArray[0, j], geometryArray[1, j], 0], lcar=meshSize))

    lines = []
    for j in range(len(points)-1):
        lines.append(geom.add_line(points[j],points[j+1]))
    lines.append(geom.add_line(points[len(points)-1],points[0]))


    #Next, we must make a closed line loop. The only argument to the add_line_loop function is a list of lines that form a closed loop. Based on the above steps, we already have that in our "lines" array and so can just pass the "lines" array to the function.
    lineLoop = geom.add_line_loop(lines)

    #Next, we define a surface using the closed line loop. The only input is the previously created line loop.
    geom.add_plane_surface(lineLoop)

    #Call the appropriate function to make the mesh. This is called on the geometry class that all of the points, lines, line loops, and surfaces have been added to.
    mesh = pygmsh.generate_mesh(geom,verbose=False)

    #path1 = os.path.join(tempfile.mkdtemp(), '.msh')
    path1 = tempfile.NamedTemporaryFile().name+'.msh'
    #path2 = os.path.join(tempfile.mkdtemp(), '.xdmf')
    path2 = tempfile.NamedTemporaryFile().name+'.xdmf'

    #meshio.write("tmp.msh", mesh)
    meshio.write(path1, mesh)

    mesh_command ="meshio-convert "+path1+" "+path2+" -p -z"
    #os.system("meshio-convert tmp.msh tmp.xdmf -p -z")

    os.system(mesh_command)

    mesh = Mesh()
    with XDMFFile(path2) as infile:
        infile.read(mesh)
    os.remove(path1)
    os.remove(path2)

    return mesh

    # Remove temporary files


class MeshModel(object):
    """
    Representation of FEM mesh for material point method using a structured
    mesh from mshr that is deformed

    Preferred method uses polygonal mesh
    """

    def __init__(self,surf_fun,bot_fun,bed_fun=None,length=4e3,Nx=10,Nz=10,xmin=0.0,dz=30.0):

        self.Nx = Nx
        self.Nz = Nz
        self.dz = dz
        self.surf_fun=surf_fun
        self.bot_fun=bot_fun
        if bed_fun==None:
                bed_fun = bot_fun
        self.bed_fun=bed_fun
        self.xmin = xmin
        self.xmax = length
        self.length = length
        self.min_grid_spacing = 0.5*length/Nx
        self.max_grid_spacing = 1.5*length/Nx
        self.nn = 501
        parameters['form_compiler']['cpp_optimize']       = True

        # Generate the mesh
        self.generate_mesh()

        # Generate the function spaces
        self.generate_function_spaces()


    def generate_mesh(self):
        """
        Generate default rectangular mesh and then deform based on surface and bottom
        topography functions
        """
        length  = self.length
        Nx = self.Nx
        Nz = self.Nz
        self.mesh = RectangleMesh(Point(0,0), Point(length, 1), Nx, Nz, "left/right")

        # Now deform top and bottom based on surface and base profiles
        coordinates = self.mesh.coordinates()
        surf = self.surf_fun(coordinates[:,0])
        bot = self.bot_fun(coordinates[:,0])
        thick = surf-bot
        coordinates[:,1] = coordinates[:,1]*thick + bot
        self.mesh.bounding_box_tree().build(self.mesh)



    def generate_function_spaces(self, order=1):
        """
        Generates the finite-element function spaces used with topological dimension :math:`d` set by variable ``self.top_dim``.
        """
        order        = 1
        space        = 'DG'
        self.Q       = FunctionSpace(self.mesh, space, order)
        self.Q2       = FunctionSpace(self.mesh, space, order+1)
        self.element = self.Q.element()
        self.Q_CG       = FunctionSpace(self.mesh, 'CG', order)
	    #We = FiniteElement("Quadrature", self.mesh.mesh.ufl_cell(), degree=order, quad_scheme='default')
        #Q = FunctionSpace(self.mesh.mesh, We)


        # topological dimension :
        self.top_dim = self.element.geometric_dimension()

        # map from verticies to nodes :
        self.dofmap  = self.Q.dofmap()

        # cell diameter :
        self.h       = project(CellDiameter(self.mesh), self.Q)  # cell diameter vector

    def get_particle_basis_functions(self, x):
        """
        Create the interpolation function associated with tracer position x needed to
        interpolate from the particles to the mesh
        """
        mesh    = self.mesh
        element = self.element

        # find the cell that contains the point (need to worry about what happens if the point is outside of
        #    the domain??)
        x_pt       = Point(*x)
        cell_id    = mesh.bounding_box_tree().compute_first_entity_collision(x_pt)
        # Check to make sure the point is in one of the cells
        if cell_id<mesh.num_cells():
            cell       = Cell(mesh, cell_id)
            coord_dofs = cell.get_vertex_coordinates()       # local coordinates

            # array for all basis functions of the cell :
            phi = np.zeros(element.space_dimension(), dtype=float)

            # compute basis function values :
            phi = element.evaluate_basis_all(x, coord_dofs, cell.orientation())


            dof = self.dofmap.cell_dofs(cell.index())
        else:
            # If the point isn't in a cell, then we set phi to zero so it doesn't count towards anything
            #   What we should do is remove the point from the array
            dof = [0,0,0]
            phi =[0.0,0.0,0.0]
        return dof, phi, cell_id#grad_phi

    def get_coords(self):
        """
        Helper function to get coordinates of the mesh
        """
        c = self.mesh.coordinates()
        return c[:,0],c[:,1]


class MeshModelPoly(MeshModel):
        """
        Representation of FEM mesh for material point method
        """

        def generate_mesh(self):
            length  = self.length
            Nx = self.Nx
            Nz = self.Nz
            thick0 = self.surf_fun(0.0)-self.bot_fun(0.0)
            dz = self.dz
            domain_vertices = list()

            # Add points along the bottom
            xptsBot = np.linspace(0,length,self.Nx)
            yptsBot = self.bot_fun(xptsBot)


            self.xptsBot = xptsBot
            self.yptsBot = yptsBot

            thick1 = self.surf_fun(length)-self.bot_fun(length)
            Nz = int(thick1/dz)
            yptsRight = np.linspace(self.bot_fun(length)+dz,self.surf_fun(length),Nz)
            xptsRight = length*np.ones(np.shape(yptsRight))

            self.xptsRight = xptsRight
            self.yptsRight = yptsRight


            xptsTop = np.linspace(length-dz,0,self.Nx-1)
            yptsTop = self.surf_fun(xptsTop)
            self.xptsTop = xptsTop
            self.yptsTop = yptsTop

            thick2 = self.surf_fun(0)-self.bot_fun(0)
            Nz = int(thick2/dz)
            yptsLeft = np.linspace(self.surf_fun(0)-dz,self.bot_fun(0)+dz,Nz)
            xptsLeft = 0.0*np.ones(np.shape(yptsLeft))

            self.xptsLeft = xptsLeft
            self.yptsLeft = yptsLeft


            xpts  =np.hstack((xptsLeft,xptsBot,xptsRight,xptsTop))
            ypts = np.hstack((yptsLeft,yptsBot,yptsRight,yptsTop))
            self.xpts = xpts
            self.ypts = ypts
            #xpts = np.append(np.append(np.append(xptsBot,xptsRight),xptsTop),xptsLeft)
            #ypts = np.append(np.append(np.append(yptsBot,yptsRight),yptsTop),yptsLeft)

            geometryArray = np.array([xpts, ypts])

            self.geometryArray = geometryArray
            # Mesh generation uses the radius so twice the mesh size
            new_mesh=meshGmsh(geometryArray, dz*2)

            #mesh = Mesh('tmp.xml')
            #mesh = Mesh()
            #with XDMFFile("tmp.xdmf") as infile:
            #    infile.read(mesh)
            self.mesh=new_mesh
            self.mesh.bounding_box_tree().build(self.mesh)

        def remesh(self,max_length = None):
            """
            Remesh cutting off mesh at max_length
            if max_length = None, then doesn't cut the domain off
            """
            dz = self.dz
            mesh = self.mesh
            bmesh = BoundaryMesh(mesh,'exterior',order=True)
            x = bmesh.coordinates()[:,0]

            if max_length == None:
                max_length = np.max(x)*10

            pts = sort_boundary_nodes(bmesh)

            # Now remove nodes that are plast the cutoff length and
            pt_new = []
            pt_flag = None
            length_flag = True
            xcliff = max_length
            for n in range(len(pts)):
                pt = pts[n]
                # We will stack x points along the calving front if they exceed the distance
                if near(pt[0],0) and pt[1]<self.bed_fun(0.0):
                        pt_new.append(pt)
                else:
                    if pt[0]<=xcliff:
                        if len(pt_new)==0:
                            pt_new.append(pt)
                        else:
                            # If there is at least one point, we calculate the distance
                            # between the new and old point
                            dist = np.sqrt((pt[0]-pt_new[-1][0])**2+(pt[1]-pt_new[-1][1])**2)
                            pt_new.append(pt)

            pt_new = np.array(pt_new)
            # The characteristic length is the radius so twice the mesh size
            new_mesh = meshGmsh(pt_new.transpose(),dz*2)


            #mesh = Mesh()
            #with XDMFFile("tmp.xdmf") as infile:
            #    infile.read(mesh)
            self.mesh=new_mesh
            self.mesh.bounding_box_tree().build(self.mesh)
            self.generate_function_spaces()
            return self
