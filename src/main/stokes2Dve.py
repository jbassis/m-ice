from fenics import *

from importlib import reload

import meshModel
reload(meshModel)
import numpy as np
from scipy.interpolate import interp1d # This is used for bed topography

import material
reload(material)
from material import *

from ufl import nabla_div

from leopart import (
    particles,
    RandomRectangle,
    l2projection,
    advect_rk3,
    advect_particles,
    assign_particle_values,
    AddDelete
)


"""
Viscoelastic implementation of quasi-material point method in fenics
"""

ffc_options = {"optimize": True, \
              "eliminate_zeros": True, \
              "precompute_basis_const": True, \
              "precompute_ip_const": True}

parameters['allow_extrapolation'] = True

def epsilon(u):
   return 0.5*(nabla_grad(u) + nabla_grad(u).T)


#def sigma(u):
#    return lambda_*nabla_div(u)*Identity(d) + 2*mu*epsilon(u)
import scipy.interpolate as interpsci
class interpND(object):
  def __init__(self, points,values):
    self.funcinterp=interpsci.LinearNDInterpolator(points,values)
    self.funcnearest=interpsci.NearestNDInterpolator(points,values)
    #self.funcRbf=Rbf(points[:,0],points[:,1],values,function='inverse')
  def __call__(self,*args):
    #t = self.funcRbf(*args)
    t=self.funcinterp(*args)
    isnan_filter = np.isnan(t)
    if np.sum(isnan_filter)==0:
      return t
    else:
      g=self.funcnearest(*args)
      t[isnan_filter]=g[isnan_filter]
    return t

def local_project(v, V, u=None):
 dv = TrialFunction(V)
 v_ = TestFunction(V)
 a_proj = inner(dv, v_)*dx
 b_proj = inner(v, v_)*dx
 solver = LocalSolver(a_proj, b_proj)
 solver.factorize()
 if u is None:
     u = Function(V)
     solver.solve_local_rhs(u)
     return u
 else:
     solver.solve_local_rhs(u)
     return

class Stokes2D:
   def __init__(self, mesh,visc_func,left_vel=0.0,right_vel=None):
       """
       Inputs
           mesh:  gridModel mesh that defines the domain and function spaces used
           visc_func:  Class that defines the viscoelastic rheology of the model
           left_vel: Horizontal velocity at the left side of the domain
           right_vel: Horizontal velocity at the right side of the domain. If set to none,
                           then a calving front/traction boundary condition is used
           elastic: Flag set to determine if the model uses elastic or a viscoelastic rheology
       """
       self.rho_i = 910.0               # kg/m3 ice density
       self.rho_w = 1020.0              # kg/m3 seawater density
       self.g = 9.81                    # m/s2
       self.f = Constant((0, -self.rho_i*self.g)) # Body Force
       self.degree = 1                  # base degree for finite elements (should be 1 or 2)

       self.t = 0.0 # Initialize time to zero


       # Use mesh to initialize geometry and a few other things
       self.mesh = mesh

       self.calving_front = False
       self.water_drag = 1e6
       self.friction = 1e8

       # Define sea level (by default it is set to zero)
       self.sea_level = 0.0

       # Regularization to enforce no-penetration condition on grounded portion of the ice sheet
       self.elastic_coeff_rock=1e6

       # Left and right horizontal velocities set at, well, left and right boundaries of the domain
       # By default the vertical velocitity is no-slip
       self.left_vel = left_vel
       self.right_vel = right_vel

       # Numerical parameters
       self.maxit = 25            # Maximum number of iterations in Picard and Newton methods
       self.tolerance = 1E-6       # Relative tolerance of solution

       # User needs to provide a callable function that computes the rheology
       self.visc_func = visc_func


       # Water pressure for portion of domain beneath sea level
       self.water_pressure = Expression(("A*(sea_level - x[1])*(x[1]<sea_level)"), A=self.rho_w * self.g,
           sea_level=self.sea_level,degree=self.degree)

       # This function determines the upward force needed to enforce the no-penetration boundary condition
       self.init_GroundPressure()

       # Boolean flag to determine if domain is beneath sea level (needed to make sure we don't apply the viscous drag force
       # when the ice is grounded)
       self.below_sea_level=Expression("x[1]<sea_level",sea_level=self.sea_level,degree=self.degree)

       # Set the mesh
       self.set_mesh(mesh)


   def effective_strain_rate_squared(self,u,Q):
      epsilon = 0.5*(nabla_grad(u) + nabla_grad(u).T)
      epsII = sqrt(epsilon[0,0]**2 + epsilon[0,1]**2)
      return epsII


   def init_GroundPressure(self):
       from scipy.optimize import fmin_cobyla
       f = self.mesh.bed_fun
       elastic_coeff_rock=self.elastic_coeff_rock
       def objective(X, P):
            x,y = X
            return np.sqrt((x - P[0])**2 + (y - P[1])**2)
       def c1(X):
            x,y = X
            return f(x) - y
       def c2(X):
            x,y = X
            return y - f(x)
       class GroundPressure(UserExpression):
            def eval(self, value, x):
                if x[1] <= f(x[0]):
                    X = fmin_cobyla(objective, x0=[x[0],f(x[0])], cons=[c1, c2], args=([x[0], x[1]],), consargs=(), disp=0)
                    value[0] = elastic_coeff_rock*objective(X, [x[0], x[1]])
                else:
                    value[0] = 0.0
            def value_shape(self):
                return ()

       self.GROUND_PRESSURE_SCALAR = GroundPressure(degree=1)

   def set_mesh(self,mesh):
       """
       This will assign the mesh to the class and initialize the function space
           Needs to be called everytime a new mesh is defined
       """
       self.mesh=mesh

       # Setup finite element function space for Taylor-Hood elements
       self.scalar_el = FiniteElement("CG", self.mesh.mesh.ufl_cell(), self.degree)    # P1 element for pressure and stress
       self.scalar_el2 = FiniteElement("CG", self.mesh.mesh.ufl_cell(), self.degree+1) # P2 element for pressure and stress
       self.vector_el = VectorElement("CG", self.mesh.mesh.ufl_cell(), self.degree+1)  # P2 element for velocity
       self.scalar = FunctionSpace(self.mesh.mesh, "CG", self.degree)                  # Function space for scalar elements
       self.vector2 = VectorFunctionSpace(self.mesh.mesh, "CG", self.degree+1)         # Needed to initialize guess for velocity
       self.vector = VectorFunctionSpace(self.mesh.mesh, "DG", self.degree)            # P1 vector function space for unit normals
       self.vector_CG = VectorFunctionSpace(self.mesh.mesh, "CG", self.degree)            # P1 vector function space for unit normals
       self.system_el = MixedElement([self.scalar_el, self.vector_el])
       self.system = FunctionSpace(self.mesh.mesh, self.system_el)                     # Pressure/Velocity/Effective Stress

   def get_coords(self):
       """
       Returns the x and z coordinates associated with the gridModel
       """
       x=self.mesh.coordinates()[:,0]
       z=self.mesh.coordinates()[:,1]
       return x,z

   def get_velocity(self):
       """
       Returns the vertext coordinates fo the velocity
       """
       ux,uz = self.u.split()
       ux = ux.compute_vertex_values()
       uz = uz.compute_vertex_values()
       return ux,uz


   def makeGroundPressure(self):
       k = self.elastic_coeff_rock
       bed = self.mesh.bed_fun
       class GroundPressure(Expression):
           def eval(self, value, x):
               dz = bed(x[0])-x[1]
               value[0] = k*dz*(dz>0)#*(dz>0)
       self.ground_pressure = GroundPressure(degree=1)

   def markBoundaries(self):
       """
       Mark the different boundaries appropriate for different boundary conditions
       The different ID correspond to the following
           1: Bottom of the glacier/shelf
           2: Right side of the domain and above the water line
           3: Left side of the domain where a velocity boundary condition is applied
           5: Right side of the domain and below the water line
       """

       sea_level = self.sea_level
       buttressing_height_min = self.buttressing_height_min
       buttressing_height_max = self.buttressing_height_max
       x,z = self.mesh.get_coords()
       left_wall = self.left_wall #np.min(x)
       right_wall = self.right_wall#np.max(x)
       top = np.max(z)
       bot = np.min(z)

       bed = self.mesh.bed_fun

       # Create MESH function over cell facets and set label to zero
       self.boundary_parts = MeshFunction('size_t', self.mesh.mesh, 1)
       self.boundary_parts.set_all(0)

       # Mark above water line as subdomain 2
       class RightTop(SubDomain):
           "Mark above water line as subdomain."
           def inside(self, x_values, on_boundary):
               "Defines boundaries of right side above water subdomain."
               return on_boundary and x_values[1] > sea_level
       GAMMA_2 = RightTop()
       GAMMA_2.mark(self.boundary_parts, 5,check_midpoint=False)


       # Mark below water boundary as subdomain 5
       class RightBelow(SubDomain):
           "Mark right side below water subdomain."
           def inside(self, x_values, on_boundary):
               "Defines boundaries of right side above water subdomain."
               return on_boundary and x_values[1] <= sea_level

       GAMMA_5 = RightBelow()
       GAMMA_5.mark(self.boundary_parts, 5,check_midpoint=False)

       # Mark below water boundary as subdomain 5
       class Buttressing(SubDomain):
           "Mark 10 meters above water line to apply buttressing"
           def inside(self, x_values, on_boundary):
               "Defines boundaries of right side above water subdomain."
               return on_boundary and (x_values[1] >= sea_level-buttressing_height_min) and (x_values[1]<=sea_level+buttressing_height_max)

       GAMMA_6 = Buttressing()
       GAMMA_6.mark(self.boundary_parts, 6,check_midpoint=False)

       # Mark bottom boundary facets as subdomain 1
       class Bottom(SubDomain):
           "Mark nodes that are in contact or below the bed"
           def inside(self, x_values, on_boundary):
               "Defines boundaries of bottom subdomain."
               return on_boundary and (x_values[1]-bed(x_values[0]) <= 10000*DOLFIN_EPS_LARGE)


       GAMMA_1 = Bottom()
       GAMMA_1.mark(self.boundary_parts, 1,check_midpoint=False)

       # Mark left velocity boundary condition
       if self.left_vel != None:
           class Left(SubDomain):
               "Mark nodes along the left wall"
               def inside(self, x_values, on_boundary):
                   "Defines boundaries of left subdomain."
                   return on_boundary and (x_values[0]-left_wall<1000000000*DOLFIN_EPS_LARGE)

           GAMMA_3 = Left()
           GAMMA_3.mark(self.boundary_parts, 3,check_midpoint=False)

       # Mark right velocity boundary condition
       if self.right_vel != None:
           class Right(SubDomain):
               "Remakes right subdomain."
               def inside(self, x_values, on_boundary):
                   "Defines boundaries of right subdomain."
                   return on_boundary and (right_wall-x_values[0]<1000000000*DOLFIN_EPS_LARGE)
           GAMMA_4 = Right()
           GAMMA_4.mark(self.boundary_parts, 4,check_midpoint=False)


       self.BCS = []
       if self.left_vel != None:
           BCL = DirichletBC(self.system.sub(1).sub(0), Constant(self.left_vel), self.boundary_parts, 3)
           #BCL = DirichletBC(self.system.sub(1), (Constant(self.left_vel),Constant(0.0)), self.boundary_parts, 3)
           self.BCS.append(BCL)

       if self.right_vel !=None:
           BCR = DirichletBC(self.system.sub(1).sub(0), Constant(self.right_vel), self.boundary_parts, 4)
           #BCR = DirichletBC(self.system.sub(1), (Constant(self.right_vel),Constant(0.0)), self.boundary_parts, 4)
           self.BCS.append(BCR)

       self.DS = Measure("ds")(subdomain_data=self.boundary_parts)


   def solve(self,p,dt=3600,tolerance=1e-6,relax_param=1.0):
       """
       Picard iteration to solve system of equations

       p is particles
       """

       # Normal vector
       z_vector = interpolate(Constant((0.0, 1.0)), self.vector2) # Vertical vector pointing upwards

       # Normal and tangent unit vectors
       N = FacetNormal(self.mesh.mesh)
       normE = as_vector([N[0], N[1]])
       tanE = as_vector([-N[1], N[0]])

       # Mark boundaries and subdomains
       self.markBoundaries()

       #Define functio space for viscosity coefficients???
       Q = self.mesh.Q

       # Extract strain and temperature as mesh functions
       # Function spaces -- should be defined once??
       Vdg = FunctionSpace(self.mesh.mesh, 'DG',1)
       Vcg = FunctionSpace(self.mesh.mesh, 'DG',1)

       # Variables to store strain and temp
       strain, temp = Function(Vdg), Function(Vcg)

       lstsq_temp = l2projection(p, Vcg, 2)
       lstsq_temp.project(temp)

       lstsq_strain = l2projection(p, Vdg, 1) # First variable???
       lstsq_strain.project_mpm(strain) # Projection is stored in phih0





       temp = self.tempModel.Temp
       u_old = Function(self.vector2)
       if self.u_k != None:
           u_k = self.u_k
           u_old.assign(u_k)
       else:
           u_k = interpolate(Constant((1E-15, 1E-15)), self.vector2)

       (q,v) = TestFunctions(self.system)
       (p, u) = TrialFunctions(self.system)

       err_rel = 1.0               # error measure ||u-u_k|| and ||p-p_k||
       count = 0                   # iteration counter
       maxit = self.maxit          # number of iterations allowed

       visc=self.visc_func



       #Q=FunctionSpace(self.mesh.mesh, "CG", 2)
       temp=interpolate(temp,Q)
       epsII=self.effective_strain_rate_squared(u_k,Q)
       eta = visc(epsII ,temp,strain,Q)
       eta_visc = self.visc_func.ductile_visc(epsII ,temp,Q)


       # Time step needs to be a Constant to avoid recompiling every time step

       water_drag = Constant(self.water_drag*self.rho_w/material.time_factor**2)*sqrt(dot(u_k,u_k))

       # Set friction coefficient
       m = self.m
       u_b = dot(tanE,u_k)
       u_b_norm = sqrt(u_b**2)
       Neff = Constant(1.0) # Don't use effective pressure in sliding law . . . for now
       tau_y = self.visc_func.cohes(strain,Q)
       friction = Neff*Constant(self.friction/material.time_factor**m)*(sqrt(u_b_norm**2 + 1e-16**2))**(m-1)
       friction = (1./friction + u_b_norm/tau_y)**(-1.0)

       lateral_drag = Constant(self.lateral_drag/material.time_factor**m)*(sqrt(u_k**2 + 1e-16**2))**(m-1)
       lateral_drag = (1./lateral_drag + sqrt(dot(u_k,u_k))/self.visc_func.yield_strength)**(-1.0)


       bcs = self.BCS
       dt_step = Constant(dt)
       elastic = Constant(self.elastic_coeff_rock)
       stokes = inner(Constant(self.rho_i/material.time_factor)*(u-u_old)/dt_step,v)*dx \
             + inner(2*eta*epsilon(u), epsilon(v))*dx \
             - inner(nabla_div(v), p)*dx \
             + inner(nabla_div(u), q)*dx \
             + inner(self.below_sea_level*self.rho_w*self.g*dot((dt_step*u)*np.abs(N[1]), z_vector)*z_vector, v)*self.DS(5) \
             + inner(elastic*dot((dt_step*u), normE), dot(v,normE))*self.DS(1) \
             + inner(self.water_pressure*normE, v)*self.DS(5) \
             + inner((self.GROUND_PRESSURE_SCALAR+self.water_pressure)*normE, v)*self.DS(1) \
             + inner(dot(v, tanE), friction*dot(u, tanE))*self.DS(1) \
             + inner(self.below_sea_level*water_drag*u, v)*self.DS(5) \
             + inner(self.water_pressure*normE, v)*self.DS(6) \
             + inner(self.below_sea_level*self.rho_w*self.g*dot((dt_step*u)*np.abs(N[1]), z_vector)*z_vector, v)*self.DS(6) \
             + inner(self.below_sea_level*water_drag*u, v)*self.DS(6) \
             + inner(v,lateral_drag*u)*dx \
             + inner(self.buttressing, v[0])*self.DS(6) \
             - inner(self.f, v)*dx \


       # Solves problem . . . .
       w = Function(self.system)
       problem = LinearVariationalProblem(lhs(stokes), rhs(stokes), w, bcs)
       solver = LinearVariationalSolver(problem)
       prm = solver.parameters
       info(prm, False)
       prm['linear_solver']            = 'mumps'
       #prm['linear_solver']            = 'lu'

       while ((float(err_rel) > tolerance) and (count < maxit)):
           solver.solve()
           self.u_p = w
           p, u = w.split(deepcopy=True)

           du = u.vector().get_local()-u_k.vector().get_local()

           err = np.linalg.norm(du)
           err_rel = err/norm(u)
           print("count = %d, relative error = %G, absolute error = %G" % (count, err_rel,err))
           alpha = self.alpha
           u.vector()[:]=alpha*u.vector()[:]+(1-alpha)*u_k.vector()[:]
           # Assign new variables to old guess
           assign(u_k, u)
           count += 1
           # Kick out if the absolute error in the velocity is less than 0.1 m/a
           if err<1e-3/material.secpera:
               break

       ux,uz = u.split()
       speed = np.sqrt(ux.compute_vertex_values()**2+uz.compute_vertex_values()**2)
       print('Max viscous speed',np.max(speed))
       epsII=self.effective_strain_rate_squared(u,Q)

       self.eta=self.visc_func(epsII ,temp,strain,Q)
       self.eta_visc=self.visc_func.ductile_visc(epsII ,temp,Q)
       self.eta_plas=self.visc_func.plastic_visc(epsII ,strain,Q)
       #self.epsII = epsII

       Q = Vdg
       epsII = Function(Q) # Plastic viscosity
       eps1 = Function(Q) # Plastic viscosity
       eps2 = Function(Q) # Plastic viscosity
       eps = epsilon(u)
       local_project(eps[0,0],Q,eps1)
       local_project(eps[0,1],Q,eps2)
       epsII.vector()[:] = np.sqrt(eps1.vector().get_local()**2+eps2.vector().get_local()**2)
       #local_project(sqrt(eps1**2 + eps2**2),Q,epsII)
       self.epsII = epsII


       if count>=maxit:
           print("WARNING: MAXIMUM NUMBER OF ITERATIONS EXCEEDED")
       self.u_k = u_k
       self.u=u
       self.u_p = w
       self.u = u
       self.p = p
       self.eta = eta

       return u,p


   def update(self,u,dt,p,remesh=True,remesh_elastic=False):
       """
       Update positions and properties of traces using 1st order Euler forward
       The RK 4 scheme is fourth order in space, but first order in time
       because we interpolate the velocities instead of solving for new
       velocities at each step
       """

       # Step 1, establish our function spaces for projections
       Q = self.mesh.Q      # Linear elements
       Q0 = FunctionSpace(self.mesh.mesh, "DG", 0) # Discontinuous elements

       epsII=self.epsII
       eta = self.eta
       eta_visc = self.eta_visc
       eta_plas = self.eta_plas

       Vdg = FunctionSpace(self.mesh.mesh, 'DG',1)


       # Variables to store strain and temp
       strain, temp = Function(Vdg), Function(Vdg)
       #lstsq_strain = l2projection(p, Vdg, 1) # First variable???
       #lstsq_strain.project_mpm(strain) # Projection is stored in phih0

       lstsq_temp = l2projection(p, Vdg, 2) # First variable???
       lstsq_temp.project(temp,253.15,273.15) # Projection is stored in phih0

       dt_min = 0.5*project(CellDiameter(self.mesh.mesh)/sqrt(dot(u, u)),Q0).compute_vertex_values()
       dt_m = np.minimum(dt,np.min(dt_min))


       #epsII = project(epsII,Vdg)
       p.interpolate(epsII,3)
       self.epsII = epsII

       (xp , pstrain , ptemp, pepsII) = (p. return_property(mesh , 0) ,
           p. return_property(mesh , 1) ,
           p. return_property(mesh , 2),
           p. return_property(mesh , 3))
       pepsII = np.maximum(pepsII,0.0)


       self.pepsII = pepsII
       self.ptemp = ptemp
       self.pstrain = pstrain
       #"""
       if self.method==1:

           deps_dt = Function(Vdg)
           deps_dt.vector()[:] = self.visc_func.strain_update(epsII.vector().get_local(),temp.vector().get_local(),strain.vector().get_local(),dt_m)
           deps_dt_eff = Function(Vdg)
           #if self.deps_dt == None:
           p.interpolate(deps_dt,1)
           pstrain_new = np.maximum(p. return_property(mesh , 1) + pstrain,0.0)
            #self.deps_dt = deps_dt

       # else:
       #     if self.deps_dt_old == None:
       #         deps_dt_old = interpolate(self.deps_dt,Vdg)
       #         deps_dt_eff.vector()[:]=1.5*deps_dt.vector().get_local()-0.5*deps_dt_old.vector().get_local()
       #         p.interpolate(deps_dt_eff,1)
       #         self.deps_dt = deps_dt
       #
       #         self.deps_dt_old = deps_dt_old
       #     else:
       #         if self.deps_dt_older==None:
       #             deps_dt_older = interpolate(self.deps_dt_old,Vdg)
       #             deps_dt_old = interpolate(self.deps_dt,Vdg)
       #             deps_dt_eff.vector()[:]=23./12*deps_dt.vector().get_local()-16./12*deps_dt_old.vector().get_local()+5./12*deps_dt_older.vector().get_local()
       #             p.interpolate(deps_dt_eff,1)
       #             self.deps_dt = deps_dt
       #             self.deps_dt_old = deps_dt_old
       #             self.deps_dt_older = deps_dt_older
       #         else:
       #             deps_dt_old = interpolate(self.deps_dt,Vdg)
       #             deps_dt_older = interpolate(self.deps_dt_old,Vdg)
       #             deps_dt_oldist = interpolate(self.deps_dt_older,Vdg)
       #             deps_dt_eff.vector()[:]=(55./24*deps_dt.vector().get_local()-59./24*deps_dt_old.vector().get_local()+37./24*deps_dt_older.vector().get_local() -9./24*deps_dt_oldist.vector().get_local())
       #             p.interpolate(deps_dt_eff,1)
       #             self.deps_dt = deps_dt
       #             self.deps_dt_old = deps_dt_old
       #             self.deps_dt_older = deps_dt_older

       else:
           pstrain_new = self.visc_func.update(pepsII,ptemp,pstrain,dt_m)
           pstrain_new = np.maximum(pstrain_new,0.0)

       #pstrain_new[xp[:,0]<1e3]=0.0
       p.change_property(pstrain_new,1)


       strain = Function(Vdg)
       lstsq_strain = l2projection(p, Vdg, 1) # First variable???
       lstsq_strain.project_mpm(strain) # Projection is stored in phih0
       #self.strain = strain


       print('Starting to remesh')
       if remesh == True:
           if remesh_elastic==False:
               #Temp = self.temp_model.advect_diffuse(T0,u,dt,self.scalar,self.boundary_parts,self.mesh.mesh)
               #self.Temp = Temp
               if self.calving_front == False:

                  # Update model mesh with new mesh
                  new_mesh=self.remesh(u,dt_m)
                  u_eff = u
                  print('Finished remeshing')
               else:
                  # Update mesh coordinates to new coordinates
                  u_elastic=self.remesh_elastic(u,dt_m)
                  #u_eff = u-u_elastic
                  #u_eff = project(u-u_elastic,self.vector2)
                  #Temp = self.temp_model.advect_diffuse(T0,u-u_elastic,dt,self.scalar,self.boundary_parts,self.mesh.mesh)
                  #self.Temp = Temp
                  ux,uz=u_elastic.split()


                  # Update mesh coordinates
                  coords = self.mesh.mesh.coordinates()
                  coords[:,0]=coords[:,0]+ux.compute_vertex_values()*dt_m
                  coords[:,1]=coords[:,1]+uz.compute_vertex_values()*dt_m

                  # And update bounding box tree
                  self.mesh.mesh.bounding_box_tree().build(self.mesh.mesh)
                  self.markBoundaries()
                  # Now remesh
                  new_mesh=self.mesh.remesh(max_length=self.mesh.length)
                  print('Finished remeshing')

               Qnew=FunctionSpace(new_mesh.mesh, "CG", 1)
               Tnew = Function(Qnew)
               self.set_mesh(new_mesh)
               length = self.mesh.length
           else:
               #xm,zm = self.tracers.get_coords()
               x,z = self.mesh.get_coords()
               xmax = np.max(x)
               u_elastic=self.remesh_elastic(u,dt_m)
               u_eff = u-u_elastic
               u_eff = project(u-u_elastic,self.vector2)
               ux,uz=u_elastic.split()


               # Update mesh coordinates
               coords = self.mesh.mesh.coordinates()
               coords[:,0]=coords[:,0]+ux.compute_vertex_values()*dt_m
               coords[:,1]=coords[:,1]+uz.compute_vertex_values()*dt_m

               # And update bounding box tree
               self.mesh.mesh.bounding_box_tree().build(self.mesh.mesh)
               self.markBoundaries()

               length = xmax*2.0
           #self.tracers.set_mesh(self.mesh)
       else:
           length = self.mesh.length


       self.tempModel.set_mesh(self.mesh.mesh)
       self.markBoundaries()
       #Temp = self.tempModel.update(u_eff,dt_m,self.boundary_parts)
       if self.u_k!=None:
           self.u_k = interpolate(self.u_k,self.vector2)




       self.strain = strain
       self.temp = temp
       self.epsII = epsII



       self.u_k = interpolate(self.u_k,self.vector2)
       #ap = advect_rk3(p, self.vector2, u, "open")
       #ap.do_step(dt_m)
       return dt_m

   def remesh_elastic(self,vel,dt):
       """
        Remesh by solving fictious elasticity problem
            We will solve an elastic problem on the mesh with deformation defined by the displacement
            associated with the velocity field along the boundaries
       Updated to use 4th order (in space) Runge-Kutta method
       """

       # Mark above water line as subdomain 2
       x,z = self.mesh.get_coords()
       left_wall = self.left_wall
       class boundary(SubDomain):
           "Mark above water line as subdomain."
           def inside(self, x_values, on_boundary):
               "Defines boundaries of right side above water subdomain."
               return on_boundary
       GAMMA_1 = boundary()
       GAMMA_1.mark(self.boundary_parts, 1)

       # Mark left boundary as subdomain 3
       class Left(SubDomain):
          "Mark nodes along the left wall"
          def inside(self, x_values, on_boundary):
              "Defines boundaries of left subdomain."
              return on_boundary and (x_values[0]-left_wall<100000000000*DOLFIN_EPS_LARGE)

       GAMMA_2 = Left()
       GAMMA_2.mark(self.boundary_parts, 3)

       # Mark left boundary as subdomain 3
       if self.right_vel !=None:
           right_wall = self.right_wall
           class Right(SubDomain):
              "Mark nodes along the right wall"
              def inside(self, x_values, on_boundary):
                  "Defines boundaries of left subdomain."
                  return on_boundary and (right_wall-x_values[0]<100000000000*DOLFIN_EPS_LARGE)
           GAMMA_3 = Right()
           GAMMA_3.mark(self.boundary_parts, 3)

       # First trial step for RK4 method
       initial_coords = self.mesh.mesh.coordinates()
       x1 = np.copy(self.mesh.mesh.coordinates()[:,0])
       z1 = np.copy(self.mesh.mesh.coordinates()[:,1])
       Q = self.vector # Function Space
       u1 = interpolate(vel, VectorFunctionSpace(self.mesh.mesh, "CG", self.degree))
       ux1,uz1=u1.split()

       # Effective velocity
       Q = VectorFunctionSpace(self.mesh.mesh, "CG", self.degree) # Function Space
       uq = interpolate(vel,Q)
       uq.vector()[:]=u1.vector().get_local()

       # Create Dirichlet boundary conditions that apply the displacement to the boundary nodes
       # Velocity boundary condition to boundary nodes that are displaced
       bc1 = DirichletBC(Q, uq, self.boundary_parts, 1)
       bc2 = DirichletBC(Q.sub(0), Constant(0.0), self.boundary_parts, 3)
       bc2.apply(uq.vector())


       new_mesh = Mesh(self.mesh.mesh)
       coords = new_mesh.coordinates()
       ux,uz = uq.split(deepcopy=True)
       coords[:,0]=coords[:,0]+ux.compute_vertex_values()*dt
       coords[:,1]=coords[:,1]+uz.compute_vertex_values()*dt

       bmesh  = BoundaryMesh(new_mesh, 'exterior')
       new_mesh = Mesh(self.mesh.mesh)
       q=ALE.move(new_mesh, bmesh)
       V = VectorFunctionSpace(new_mesh, 'CG', 1)
       q = interpolate(q,V)
       q.vector()[:] /= dt
       self.q = q
       return q

   def remesh(self,u,dt):
       Q = self.mesh.Q_CG
       ux,uz = u.split()
       ux = project(ux,Q).vector().get_local()
       uz = project(uz,Q).vector().get_local()


       # To do this we create a function and then mark top,bot, left and right
       boundary_dofs= Function(Q)
       bcTop = DirichletBC(Q, 12, self.boundary_parts, 2)
       bcBot = DirichletBC(Q, 14, self.boundary_parts, 5)
       bcBed = DirichletBC(Q, 14, self.boundary_parts, 1)
       bcTop.apply(boundary_dofs.vector())
       bcBot.apply(boundary_dofs.vector())
       bcBed.apply(boundary_dofs.vector())

       # Extract coordinates corresponding to dofs
       dof_coords = Q.tabulate_dof_coordinates().reshape((-1, 2));xdof=dof_coords[:,0];zdof=dof_coords[:,1]

       # Save dof coordinates for debuggins and plotting
       self.xdof = xdof
       self.zdof = zdof

       # Extract coordinates of bottom of mesh
       bot_nodes = boundary_dofs.vector().get_local()==14
       xbot=xdof[bot_nodes]
       zbot=zdof[bot_nodes]

       # Extract velocity components at bottom of the mesh
       ux_bot = ux[bot_nodes]
       uz_bot = uz[bot_nodes]

       # Update bottom nodes
       xbot+=ux_bot*dt
       zbot+=uz_bot*dt

       # Now update top nodes
       # Extract coordinates of top of mesh
       top_nodes = boundary_dofs.vector().get_local()==12
       xtop=xdof[top_nodes]
       ztop=zdof[top_nodes]

       # Extract velocity components at top of the mesh
       ux_top = ux[top_nodes]
       uz_top = uz[top_nodes]

       # Update top nodes
       xtop+=ux_top*dt
       ztop+=uz_top*dt

       # Create interpolation function to define top and bottom of the ice
       id1 = np.argsort(xbot)
       id2 = np.argsort(xtop)
       bot_fun= interp1d(xbot[id1],zbot[id1],fill_value="extrapolate")
       surf_fun = interp1d(xtop[id2],ztop[id2],fill_value="extrapolate")

       self.xtop = xtop[id2]
       self.ztop = ztop[id2]
       self.xbot = xbot[id1]
       self.zbot = zbot[id1]


       # Create a new mesh
       bed_fun = self.mesh.bed_fun
       Nx = self.mesh.Nx
       Nz = self.mesh.Nz
       length = self.mesh.length
       new_mesh = gridModel.MeshModelPoly(surf_fun,bot_fun,bed_fun,Nx=Nx,Nz=Nz,length=length)

       return new_mesh
