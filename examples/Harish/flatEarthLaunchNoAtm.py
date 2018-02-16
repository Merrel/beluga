
import numpy as np

import beluga.Beluga as Beluga
import beluga.bvpsol as bvpsol
import beluga.bvpsol.algorithms as algorithms
import beluga.optim.Problem
from beluga.optim.problem import *
from beluga.continuation import *
import logging

def get_problem():
    """Flat Earth launch example with no atmosphere."""

    # Rename this and/or move to optim package?
    problem = beluga.optim.Problem('flatEarthLaunchNoAtm')

    # Switch off DAE mode
    problem.mode = 'analytical'

    # Define independent variables
    problem.independent('t', 's')

    # Define equations of motion
    problem.state('x', 'vx', 'm')   \
           .state('y', 'vy', 'm')   \
           .state('vx', 'T*cos(pitch)/M', 'm/s') \
           .state('vy', '(T*sin(pitch)/M) - g', 'm/s')
    # Define controls
    problem.control('pitch','rad')

    # Define costs
    problem.cost['path'] = Expression('1','s')

    # Define constraints
    problem.constraints('default',0) \
                        .initial('x-x_0','m')    \
                        .initial('y-y_0','m')    \
                        .initial('vx-vx_0','m/s')  \
                        .initial('vy-vy_0','m/s') \
                        .terminal('y-y_f','m')   \
                        .terminal('vx-vx_f','m/s')
#                        .terminal('vy-vy_f','m/s')
    # problem.constraints().interior_point('(x-x1)^2+(y-y1)^2','m^2')

    # Define constants (change to have units as well)
    problem.constant('g','9.81','m/s^2') \
           .constant('M','100000','kg') \
           .constant('T','500000','n')
    # problem.constant('x1','7.5','m')
    # problem.constant('y1','-15','m')

    # problem.quantity('gDown','g*sin(theta)')

    problem.scale.unit('m',1)     \
                   .unit('s',1)\
                   .unit('kg',1)   \
                   .unit('rad',1) \
                   .unit('n',1)

    problem.bvp_solver = algorithms.MultipleShooting(derivative_method='fd',tolerance=1e-4, max_iterations=1000, verbose = True, cached=False, number_arcs=4, max_error=100)
    #problem.bvp_solver = algorithms.SingleShooting(derivative_method='fd',tolerance=1e-4, max_iterations=50, verbose = True, cached=False)
    # problem.bvp_solver = algorithms.BroydenShooting(tolerance=1e-4, max_iterations=1000)

    # Can be array or function handle
    # TODO: implement an "initial guess" class subclassing Solution
    # problem.guess = bvpsol.bvpinit(np.linspace(0,1,2), [0,0,1,-0.1,-0.1,-0.1,0.1])
    # problem.guess.parameters = np.array([0.1,0.1,0.1,0.1,0.1])
    problem.guess.setup('auto',
                    start=[0,0,0,0],          # Starting values for states in order
                    direction='forward',
                    costate_guess = 0.1#[-0.1,-0.1,-0.1,-0.1]
                    )

    # Figure out nicer way of representing this. Done?
    problem.steps.add_step('bisection') \
                    .num_cases(1001) \
                    .terminal('y', 150000) \
                    .terminal('vx', 7800) \
                    .terminal('vy', 0)

    # (
    # problem.steps.add_step().num_cases(2)
    #                  .terminal('x', 30.0)
    #                  .terminal('y',-30.0),

    # problem.steps.add_step()
    #                 .num_cases(10)
    #                 .terminal('x', 1000.0)
    #                 .terminal('y',-1000.0)
    # )
    return problem

if __name__ == '__main__':
    Beluga.run(get_problem(), display_level=logging.DEBUG)
