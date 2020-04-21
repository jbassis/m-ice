"""
Set of routines to solve temperature advection and diffusion to avoid keeping temperature as
a tracer particle
"""
from fenics import *
import numpy as np
from material import secpera
from scipy.special import erfc
import pickle
parameters['allow_extrapolation'] = True

advect_code = '''
#include <pybind11/pybind11.h>
#include <math.h>
class advect : public UserExpression
{
public:

void eval(Array<double>& values, const Array<double>& x) const {
    // Evaluate velocity at x and store in array u
    Array<double> u(2);
    velocity->eval(u,x);


    // Use method of characteristics to find x_old = x-u*dt
    Array<double> P_old(2);
    P_old[0] = x[0]-u[0]*dt;
    P_old[1] = x[1]-u[1]*dt;

    // cell_id is negative if point is outside of cell
    int cell_id = tree->compute_first_collision(P_old);
    if (cell_id<0) {
        temp->eval(values, x);
    }
    else{
        // Evaluate strain at upstream grid point
        temp->eval(values, P_old);

    }
}

std::shared_ptr<const Function> temp;
std::shared_ptr<const Function> velocity;
std::shared_ptr<const BoundingBoxTree> tree;
double dt;
};
'''


class tempModel:
    def __init__(self,mesh=None,Tb=0,Ts=-20,beta=0.0025):
        # Diffusion coefficient
        self.kappa = 1.1788e-6*secpera
        # Surface and bottom temperatures in Kelvin
        self.Tb = 273.15+Tb
        self.Ts = 273.15+Ts

        # Define function space parameters
        self.degree = 1   #Degree
        self.space = 'CG' #Type of element
        if mesh !=None:
            self.set_mesh(mesh)

        # Define parameters for temperature distribution to be used as a
        #  boundary condition on the left (and/or right)
        #  side of the domain
        self.beta = 0.0025 # Transition to bed temp


    def set_mesh(self,mesh):
        self.mesh = mesh
        self.Q = FunctionSpace(mesh, self.space, self.degree)


    def set_temp(self,x,surf_fun,bot_fun,bed_fun):
        """
        Define initial and incoming temperature field where
         xi = (z-zb)/H and;
         Tinit = erf(xi/sqrt(2*beta))*(Tb-Ts)+Ts
         beta is a parameter that controls how quickly the temperature
         transitions to Tb
        """
        beta = self.beta
        Tb=self.Tb
        Ts=self.Ts
        mesh = self.mesh

        gdim = self.mesh.geometry().dim()
        Q = self.Q
        x = Q.tabulate_dof_coordinates().reshape((-1, gdim))
        Tinit = Function(Q)
        Temp = Function(Q)
        zb = bot_fun(x[:,0])
        thick = surf_fun(x[:,0])-bot_fun(x[:,0])
        xi = (x[:,1]-zb)/thick
        temp_vals = erfc(xi/sqrt(2*beta))*(Tb-Ts)+Ts
        #temp_vals = Tb + xi(Ts-Tb)
        Tinit.vector().set_local(temp_vals)
        Temp.vector().set_local(temp_vals)

        self.Temp=Temp
        self.incoming_temp=Tinit

    def load_temp(self,file_name):
            self.Temp = Function(self.Q,file_name)

    def save(self,fname):
        param_dict = {'kappa':self.kappa,'Tb':self.Tb,'Ts':self.Ts,
                        'degree':self.degree,'beta':self.beta,'space':self.space}
        f = open(fname,"wb")
        pickle.dump(param_dict,f)
        f.close()


    def load(self,fname_params):
        # Load pickled parameters
        f = open(fname_params,"rb")
        params = pickle.load(f)
        f.close()
        for key in params:
            setattr(self, key, params[key])


    def advect(self,u,dt):
        """
        Using C++ Expression to advect temperature
        Input: u: velocity field
               dt:time step
        """
        f = Expression(advect_code, element=FiniteElement(self.space, self.mesh.ufl_cell(), self.degree))
        # We pass the tree as an argument so that we can tell
        #   if the point is in the mesh or not
        self.mesh.bounding_box_tree().build(self.mesh)
        tree = BoundingBoxTree()
        tree.build(self.mesh)
        f.tree=tree
        # Initial value of temperature
        f.temp = self.Temp
        # Time step
        f.dt = dt
        # Velocity
        f.velocity = u
        # Create a function in function space Q (this is where our dofs live)
        T_new = Function(self.Q)
        # Interpolate temperature to the new function space
        #  it is in the interpolation step where we introduce
        #  errors, but these errors should be small
        T_new.interpolate(f)
        self.Temp = T_new
        return T_new



    def update(self,velocity,dt,boundary_parts):
        # Take one advection step
        #Temp = self.advect(velocity,dt)
        Temp = self.Temp

        # Apply boundary conditions to advected temperature to
        #   set temp at edges to initial values
        BC_left = DirichletBC(self.Q, self.incoming_temp, boundary_parts, 3)
        BC_bottom = DirichletBC(self.Q, Constant(self.Tb), boundary_parts, 1)
        bcs = [BC_bottom,BC_left]
        for bc in bcs:
            bc.apply(Temp.vector())

        # In case there are any advected temperatures that are incorrect . . .
        gdim = self.mesh.geometry().dim()
        Q = self.Q
        x_dof_temp = Q.tabulate_dof_coordinates().reshape((-1, gdim))
        idz=Temp.vector().get_local()==0
        nn = sum(idz)
        if nn>0:
             xz=self.mesh.coordinates()
             values = Temp.vector().get_local()
             not_zero = values>0.0
             interp_fun = interpND(x_dof_temp[not_zero],values[not_zero])
             interp_vals = interp_fun(x_dof_temp[idz,0],x_dof_temp[idz,1])
             values[idz]=interp_vals
             Temp.vector().set_local(values)

        temp_vals = Temp.vector().get_local()
        temp_vals = np.minimum(temp_vals,self.Tb)
        temp_vals = np.maximum(temp_vals,self.Ts)
        Temp.vector().set_local(temp_vals)
        # Now diffuse and apply boundary conditions to temperature
        Temp = self.diffuse(dt,boundary_parts)
        temp_vals = Temp.vector().get_local()
        temp_vals = np.minimum(temp_vals,self.Tb)
        temp_vals = np.maximum(temp_vals,self.Ts)
        Temp.vector().set_local(temp_vals)
        self.Temp = Temp
        return Temp


    def diffuse(self,dt,boundary_parts):
        """
        Function updates diffusion using an implicit time step
        # T = T + dt*kappa*(T_xx+T_zz)
        Inputs:
            Tinit: Initial temperature field
            dt : time step
            Q : Function space
            boundary_parts : labeled portions of the boundary
        Output:
            T : Temperature at the end of the time step
        """
        T0 = self.Temp
        mesh = self.mesh
        Q = self.Q
        T = Function(Q)

        # Test and trial functions
        phi, v = TrialFunction(Q), TestFunction(Q)

        BC_left = DirichletBC(Q, self.incoming_temp, boundary_parts, 3)
        BC_bottom = DirichletBC(Q, Constant(self.Tb), boundary_parts, 1)
        bcs = [BC_bottom,BC_left]

        # Bilinear form
        F = inner(phi,v)*dx -inner(T0,v)*dx + Constant(0.5*dt*self.kappa)*inner(grad(phi),grad(v))*dx \
            + Constant(0.5*dt*self.kappa)*inner(grad(T0),grad(v))*dx
        #a = lhs(F)
        #L = rhs(F)

        #A = assemble(a)
        #b = assemble(L)
        #for bc in bcs:
        #    bc.apply(A,b)
        #solver = solver(A)
        #solver.solve(T.vector(), b)
        problem = LinearVariationalProblem(lhs(F),rhs(F), T, bcs)
        solver = LinearVariationalSolver(problem)
        solver.solve()
        #T.assign(T+dT)
        return T




    def update_SUGG(self,velocity,dt,boundary_parts):
        """
        SUPG stabilized advection-diffusion
        """
        mesh = self.mesh
        T0 = self.Temp

        # Function space that we work in
        Q = self.Q

        # Define function for temperature
        T = Function(Q)

        # Get cell size for stabilization
        h = CellSize(mesh)

        # Test and trial functions
        u, v = TrialFunction(Q), TestFunction(Q)

        # Mid-point solution
        u_mid = 0.5*(T0 + u)

        # Residual
        r = u - T0 + dt*(dot(velocity, grad(u_mid)) - self.kappa*div(grad(u_mid)))

        # Galerkin variational problem
        F = v*(u - T0)*dx + dt*(v*dot(velocity, grad(u_mid))*dx \
                              + self.kappa*dot(grad(v), grad(u_mid))*dx)

        # Add SUPG stabilisation terms
        vnorm = sqrt(dot(velocity, velocity))
        Pe = 1.0/self.kappa
        penalty = 0.5*h*pow(4.0/(Pe*h)+2.0*vnorm,-1.0)

        F += penalty*dot(velocity, grad(v))*r*dx


        # Create bilinear and linear forms
        a = lhs(F)
        L = rhs(F)

        # Set incoming temperature to the incoming temperature profile
        BC_left = DirichletBC(Q, self.incoming_temp, boundary_parts, 3)

        # Set bottom temperature to Tb (presumably the pressure melting temp)
        BC_bottom = DirichletBC(Q, Constant(self.Tb), boundary_parts, 1)

        bcs = [BC_left,BC_bottom]

        # Assemble and solve the system
        A = assemble(a)
        b = assemble(L)
        for bc in bcs:
            bc.apply(A,b)
        solver = LUSolver(A)
        solver.solve(T.vector(), b)
        temp_vals = T.vector().get_local()
        temp_vals = np.minimum(temp_vals,self.Tb)
        temp_vals = np.maximum(temp_vals,self.Ts)
        T.vector().set_local(temp_vals)
        self.Temp = T
        return T

class temp_init(UserExpression):
    def __init__(self, Ts, Tb,surf_fun, bed_fun, **kwargs):
        super().__init__(**kwargs)
        self.Ts = Ts
        self.Tb = Tb
        self.surf= surf_fun
        self.bed = bed_fun
        self.beta = 0.0025

    def eval(self, value,x):
        xi = (x[1]-self.bed(x[0]))/(self.surf(x[0])-self.bed(x[0]))
        value[0] = erfc(xi/sqrt(2*self.beta))*(self.Tb-self.Ts)+self.Ts+ 273.15
        return value

    def value_shape(self):
        return ()

strain_init = Expression("0.0", degree=1)
