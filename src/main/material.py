
import numpy as np
from fenics import *

time_factor = 86400.0*365.24
secpera = 86400.0*365.24

class glenFlow(object):
    #Computes harmonic mean of diffusion and dislocation creep based ice rheology
    R = 8.314
    def __init__(self,grain_size=1e-3,*args,**kwargs):
        # Parameters in Glen's Flow Law Rheology (from Cuffey and Paterson)
        self.enhancement = 1.0
        self.E1 = 60e3  # Activation of ice warmer than -10
        self.E2 = 115e3 # Activation energy of ice colder than -10
        self.A  = 3.5e-25 # Prefactor
        self.pre1 = self.enhancement*(np.exp(self.E1/self.R/263.15)*self.A)**(-1.0/3)/(time_factor)**(1./3)
        self.pre2 = self.enhancement*(np.exp(self.E2/self.R/263.15)*self.A)**(-1.0/3)/(time_factor)**(1./3)
        self.set_grain_size(grain_size)
        self.E_diff = 59.4e3

        # Additional parameterself.B0 = 697.3/(time_factor)**(1./3)
        self.plastic = True
        self.yield_strength = 200e3
        self.yield_min = 10e3
        self.crit_strain = 0.1
        self.visc_min = 1e8
        self.visc_max = 1e18
        self.mu = 0.0
        self.num_yielded = 0.0

    def set_grain_size(self, grain_size):
        self.grain_size = grain_size
        self.pre_diff = (1.2e-10/(grain_size**2))**(-1.0)/(time_factor)

    def ductile_visc(self,epsII,tempFun,functionSpace=None):
        """
        Calculate effective viscosity
        Input: u_values: velocity fields (as Fenics array thingy)
               scalar : the appropriate FENICS function space we live in
        Return: scalar effective viscosity in units of Pa s
        """

        # If no function space is provided then these are just raw nparrays
        if functionSpace == None:
            temp = tempFun
            epsII = np.maximum(epsII,0.0)
            Bdiff =  self.pre_diff*np.exp(self.E_diff/(8.314*temp))
            Bdisl = (self.pre1*np.exp(self.E1/(3.0*8.314*temp))*(temp<=263.15) + self.pre2*np.exp(self.E2/(3.0*8.314*temp))*(temp>263.15))
            visc = 0.5*(epsII**(2.0/3.0)/Bdisl + 1.0/Bdiff)**(-1.0)
        else:
            Bdisl = Function(functionSpace)
            Bdiff = Function(functionSpace)
            temp = interpolate(tempFun,functionSpace).vector().get_local()
            Bdiff.vector()[:] =  self.pre_diff*np.exp(self.E_diff/(8.314*temp))
            Bdisl.vector()[:] = (self.pre1*np.exp(self.E1/(3.0*8.314*temp))*(temp<=263.15) + self.pre2*np.exp(self.E2/(3.0*8.314*temp))*(temp>263.15))
            visc = 0.5*(epsII**(2.0/3.0)/Bdisl + 1.0/Bdiff)**(-1.0)
        return visc

    def cohes(self,strainFun,functionSpace=None):
        if functionSpace == None:
            strain = strainFun
            tau_y = np.maximum(self.yield_strength -(self.yield_strength -self.yield_min)*strain/self.crit_strain,self.yield_min)
            assert(np.min(tau_y)>0)
        else:
            tau_y = Function(functionSpace)
            strain = interpolate(strainFun,functionSpace).vector().get_local()
            tau_y.vector()[:] = np.maximum(self.yield_strength -(self.yield_strength -self.yield_min)*strain/self.crit_strain,self.yield_min)
            assert(np.min(tau_y.vector().get_local())>0)
        return tau_y

    def plastic_visc(self,epsII,strain,functionSpace=None):
        """
        Calculate plastic effective viscosity
        Input:
            epsII_squared: second invariant squared in UFL form
            strain: plastic strain on markers as an interpolation object or expression as a FENICS function
            functionSpace: FENICS function space for cohesion
            epsII_visc: Viscous part of epsII in UFL form

        Return: scalar effective viscosity in units of Pa s in FENICS UFL form
        """

        tau_y = self.cohes(strain,functionSpace)
        if functionSpace ==None:
            visc = 0.5*tau_y/(epsII+(1e-16))
        else:
            visc = 0.5*tau_y/(epsII+Constant(1e-16))
        return visc

    def update(self,epsII,temp,strain,dt):
        """ Needs to be called at the particle level """

        eta_visc = self.ductile_visc(epsII,temp,functionSpace=None)
        eta_plas = self.plastic_visc(epsII,strain,functionSpace=None)
        strain_new = strain + (eta_visc>eta_plas)*np.maximum(epsII,0.0)*dt
        return strain_new

    def strain_update(self,epsII,temp,strain,dt):
        """ Needs to be called at the particle level """

        eta_visc = self.ductile_visc(epsII,temp,functionSpace=None)
        eta_plas = self.plastic_visc(epsII,strain,functionSpace=None)
        #eta = (1/eta_visc + 1/eta_plas)**(-1)
        #eta = self.__call__(epsII,temp,strain,functionSpace=None)
        #F = 2*eta*epsII
        #epsII_visc = F/eta_visc
        epsII_visc = 0.0 # In case we want to subtract viscous strain rate
        deps =  (eta_visc>eta_plas)*np.maximum(epsII,0.0)*dt
        self.num_yielded = sum(eta_visc>eta_plas)/len(epsII)
        return deps

    def __call__(self,epsII,temp,strain,functionSpace=None):
        """
        Calculate effective viscosity
        Input: u_values: velocity fields (as Fenics array thingy)
               scalar : the appropriate FENICS function space we live in
        Return: scalar effective viscosity
        """
        visc = self.ductile_visc(epsII,temp,functionSpace)
        if self.plastic==True:
            tau_y = self.cohes(strain,functionSpace)
            visc =  Constant(0.5*self.visc_min/time_factor) + (1.0/visc + 2*epsII/tau_y + 2.0/Constant(self.visc_max/time_factor))**(-1.0)
        return visc
