import numpy as np
import beluga.bvpsol as bvpsol
import beluga.bvpsol.algorithms as algorithms
import beluga.optim.Problem
from beluga.optim.problem import *
from beluga.continuation import *
from math import *

import functools

def get_problem():
    # Figure out way to implement caching automatically
    #@functools.lru_cache(maxsize=None)


    #
    # from beluga.utils import keyboard
    # keyboard()

    """A simple planar hypersonic problem example."""

    # Rename this and/or move to optim package?
    problem = beluga.optim.Problem('planarHypersonic')
    problem.mode='analytical'
    # problem = beluga.optim.Problem()

    # Define multi-phase problem
    # problem.numPhases = 1


    # Define independent variables
    problem.independent('t', 's')

    # Define quantities used in the problem
    problem.quantity('rho','rho0*exp(-h/H)')
    problem.quantity('Cd','1.2')
    problem.quantity('D','0.5*rho*v^2*Cd*Aref')
    problem.quantity('r','re+h')

    # Define equations of motion
    problem.state('h','v*sin(gam)','m')   \
           .state('theta','v*cos(gam)/r','rad')  \
           .state('v','-D/mass - mu*sin(gam)/r**2','m/s') \
           .state('gam','(v/r - mu/(v*r^2))*cos(gam)','rad')

    # Define controls
    problem.control('alpha','rad')

    # Define costs
    problem.cost['path'] = Expression('(alpha^2)+((1)*(v**3))','m^3/s^3')

    # Define constraints
    problem.constraints().initial('h-h_0','m') \
                        .initial('theta-theta_0','rad') \
                        .initial('v-v_0','m/s') \
                        .terminal('h-h_f','m')  \
                        .terminal('theta-theta_f','rad')

    # Define constants
    problem.constant('mu', 398600.4418*1e9, 'm^3/s^2') # Gravitational parameter, m^3/s^2
    problem.constant('rho0', 1.225, 'kg/m^3') # Sea-level atmospheric density, kg/m^3
    problem.constant('H', 8500, 'm') # Scale height for atmosphere of Earth, m

    problem.constant('mass',2000,'kg') # Mass of vehicle, kg
    problem.constant('re',6378137,'m') # Radius of planet, m
    problem.constant('Aref',pi*(1.25)**2,'m^2') # Reference area of vehicle, m^2

    problem.bvp_solver = algorithms.MultipleShooting(derivative_method='fd',tolerance=1e-4, max_iterations=1000, verbose = True, cached=False, number_arcs=4, max_error=100)
    #problem.bvp_solver = algorithms.SingleShooting(derivative_method='fd',tolerance=1e-4, max_iterations=1000, verbose = True, cached = False, number_arcs=2, max_error=200)
    #problem.bvp_solver = algorithms.SingleShooting(derivative_method='fd',tolerance=1e-4, max_iterations=100000, verbose = True, cached = False)

    problem.scale.unit('m','h')         \
                   .unit('s','h/v')     \
                   .unit('kg','mass')   \
                   .unit('rad',1)

    problem.guess.setup('auto',start=[120000,0,7800,-90*pi/180],costate_guess=-0.1)
    #problem.guess.setup('auto',start=[80000,3.38575809e-21,5000,7.98617365e-02],direction='forward',time_integrate=229.865209,costate_guess =[-1.37514494e+01,3.80852584e+06,-3.26290152e+03,-2.31984720e-14])
    # Figure out nicer way of representing this. Done?

    problem.steps.add_step().num_cases(101) \
                            .terminal('h', 0)
 #                           .terminal('theta', 500/6378.137)

    problem.steps.add_step().num_cases(31)  \
                            .terminal('theta', 500/6378.137)
    # #
    # problem.steps.add_step()
    #                 .num_cases(3)
    #                 .terminal('x', 40.0)
    #                 .terminal('y',-40.0)
    # )
    return problem

if __name__ == '__main__':
    import beluga.Beluga as Beluga
    problem = get_problem()
    sol = Beluga.run(problem)
