"""
Remove self intersecting nodes manually
"""
from   fenics import *
import numpy  as np
from mshr import *
import os
import math
import pygmsh
#Meshio for saving mesh to file for use in fenics
import meshio

import mesh_functions
reload(mesh_functions)
from mesh_functions import *


dz = model.mesh.dz
bmesh = BoundaryMesh(model.mesh.mesh,'exterior',order=True)
pt = sort_boundary_nodes(bmesh)


# Remove two points that must be determined by manual inspection
pts = vstack((pt[0:461,:],pt[465::,:]))## Identify which points to remove
#pts = vstack((pt[0:483,:],pt[294:394,:],pt[397::,:]))## Identify which points to remove

#pts = vstack((pts[0:369,:],pts[414::,:]))## Identify which points to remove


# Make a new mesh based on remeshing stuff . . .
# Now remove nodes that are plast the cutoff length and
pt_new = []
pt_flag = None
length_flag = True
xcliff = max_length
for n in range(len(pts)):
    pt = pts[n]
    # We will stack x points along the calving front if they exceed the distance
    if near(pt[0],0) and pt[1]<model.mesh.bed_fun(0.0):
        #if len(pt_new)==0:
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
meshGmsh(pt_new.transpose(),dz*2)


mesh = Mesh()
with XDMFFile("tmp.xdmf") as infile:
    infile.read(mesh)

model.mesh.mesh=mesh
model.mesh.mesh.bounding_box_tree().build(model.mesh.mesh)
model.mesh.generate_function_spaces()

model.set_mesh(model.mesh)
model.tracers.set_mesh(model.mesh)
model.tempModel.set_mesh(model.mesh.mesh)
model.tracers.update_tracer_interp_functions()
print('Updated tracer interpolation function')
node_vars = model.tracers.tracers_to_nodes()
print('Starting to remove and reseed')
model.tracers.remove_and_reseed(node_vars,0,length,model.mesh.surf_fun,model.mesh.bed_fun)
print('Updating tracer interpolation functions')
model.tracers.update_tracer_interp_functions()
print('Finished updating tracer interpolation functions')
#model.u_k = None
model.u_k = interpolate(model.u_k,model.vector2)

# Now update everything else
