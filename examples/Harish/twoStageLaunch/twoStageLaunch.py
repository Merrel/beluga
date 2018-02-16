
from beluga.bvpsol.algorithms.rashs import *
from math import *

multiPhaseProblem = rashs('twoStageLaunch')

# Define states
multiPhaseProblem.state('h','m')   \
       .state('theta','rad')  \
       .state('v','m/s') \
       .state('gam','rad')

# Define the equations of motion for each phase
multiPhaseProblem.defineDynamicsFunction(['v*sin(gam)','v*sin(gam)'])
multiPhaseProblem.defineDynamicsFunction(['v*cos(gam)/r','v*cos(gam)/r'])
multiPhaseProblem.defineDynamicsFunction(['(T1*cos(alpha)/mass1)-D/mass1 - mu*sin(gam)/r**2','(T2*cos(alpha)/mass2)-D/mass2 - mu*sin(gam)/r**2'])
multiPhaseProblem.defineDynamicsFunction(['(T1*sin(alpha)/(mass1*v))+((v/r - mu/(v*r^2))*cos(gam))','(T2*sin(alpha)/(mass2*v))+((v/r - mu/(v*r^2))*cos(gam))']) # phase 1

# Define costs (units have to be same for all)
# multiPhaseProblem.pathCost('(alpha^2)+(0.0001*sqrt(rho)*(v**3))','m/s^3')
# multiPhaseProblem.pathCost('(alpha^2)+(0.0001*sqrt(rho)*(v**3))','m/s^3')

multiPhaseProblem.pathCost('1','nd')
multiPhaseProblem.pathCost('1','nd')

# Define scalar switching condition for each phase [g1,g2,g3,...,gn]
multiPhaseProblem.switchingConditionSetup(['(t-ts)']) # phase 1
multiPhaseProblem.switchingConditionSetup(['(-t+ts)']) # phase 2

multiPhaseProblem.rashsDynamicsSetup()
multiPhaseProblem.rashsPathCostSetup()

ocp = beluga.OCP(multiPhaseProblem.name)

# Define independent variables
ocp.independent('t', 's')

# Define quantities used in the problem
ocp.quantity('rho','rho0*exp(-h/H)')
ocp.quantity('D','0.5*rho*v^2*Cd*Aref')
ocp.quantity('r','re+h')

# Define equations of motion
ocp.state(multiPhaseProblem.stateVariable[0],multiPhaseProblem.rashsDynamics[0],multiPhaseProblem.stateUnits[0])   \
   .state(multiPhaseProblem.stateVariable[1],multiPhaseProblem.rashsDynamics[1],multiPhaseProblem.stateUnits[1])  \
   .state(multiPhaseProblem.stateVariable[2],multiPhaseProblem.rashsDynamics[2],multiPhaseProblem.stateUnits[2]) \
   .state(multiPhaseProblem.stateVariable[3],multiPhaseProblem.rashsDynamics[3],multiPhaseProblem.stateUnits[3])
print(multiPhaseProblem.rashsDynamics[1])
# Define controls
ocp.control('alpha','rad')

# Define costs
# print(multiPhaseProblem.rashsPathCost)
ocp.path_cost(multiPhaseProblem.rashsPathCost[0],multiPhaseProblem.pathCostUnits[0])
#problem.cost['path'] = Expression('(alpha^2)+((1)*(v**3))','m^3/s^3')


# Define constraints
ocp.constraints().initial('h-h_0','m') \
                 .initial('theta-theta_0','rad') \
                 .initial('v-v_0','m/s') \
                 .terminal('h-h_f','m')  \
                 .terminal('v-v_f','m/s') \
                 .terminal('gam-gam_f','rad')

ocp.scale(m='h',s='h/v',kg='mass1',rad=1,nd=1)

# Define constants
ocp.constant('mu', 398600.4418e9, 'm^3/s^2') # Gravitational parameter, m^3/s^2
ocp.constant('rho0', 1.125, 'kg/m^3') # Sea-level atmospheric density, kg/m^3
ocp.constant('H', 8500.0, 'm') # Scale height for atmosphere of Earth, m
ocp.constant('mass1',50000.0,'kg') # Mass of vehicle, kg
ocp.constant('mass2',10000.0,'kg') # Mass of vehicle, kg
ocp.constant('re',6378137.0,'m') # Radius of planet, m
ocp.constant('Aref',pi*(1.25)**2,'m^2') # Reference area of vehicle, m^2
ocp.constant('slpRashs',0.001,'nd')
ocp.constant('Cd',0.2,'nd')
ocp.constant('T1',1500000,'kg*m/s^2')
ocp.constant('T2',100000,'kg*m/s^2')
ocp.constant('ts',250,'s')

bvp_solver = beluga.bvp_algorithm('MultipleShooting',
                        derivative_method='fd',
                        tolerance=1e-4,
                        max_iterations=400,
                        verbose = True,
                        max_error=10000,
                        # number_arcs=2
             )

#problem.bvp_solver = algorithms.MultipleShooting(derivative_method='fd',tolerance=1e-4, max_iterations=1000, verbose = True, cached=False, number_arcs=4, max_error=100000)
#problem.bvp_solver = algorithms.SingleShooting(derivative_method='fd',tolerance=1e-4, max_iterations=1000, verbose = True, cached = False, number_arcs=2, max_error=200)
#problem.bvp_solver = algorithms.SingleShooting(derivative_method='fd',tolerance=1e-4, max_iterations=100000, verbose = True, cached = False)

guess_maker = beluga.guess_generator('auto',
                start=[0.1,0.1,1,40*pi/180],          # Starting values for states in order
                direction='forward',
                costate_guess = [-0.1,0.1,0.1,0]
)

continuation_steps = beluga.init_continuation()

# problem.guess.setup('auto',start=[120000.0,0,6000.0,-0*pi/180],costate_guess=-0.1)
# problem.guess.setup('auto',start=[120000.0,0,6000.0,-0*pi/180],direction='forward',time_integrate=0.1,costate_guess =[-0.1,0.1,0.1,0])
# Figure out nicer way of representing this. Done?

continuation_steps.add_step('bisection') \
                .num_cases(200) \
                .terminal('h', 150000) \
                .terminal('v', 7800) \
                .terminal('gam', 0)
# continuation_steps.add_step('bisection') \
#                 .num_cases(200) \
#                 .terminal('theta', 350/3396.2)

# continuation_steps.add_step('bisection') \
#                 .num_cases(1000) \
#                 .const('Cd2', 14)

# continuation_steps.add_step('bisection') \
#                 .num_cases(1000) \
#                 .const('Cd2', 14)

beluga.solve(ocp,
             method='traditional',
             bvp_algorithm=bvp_solver,
             steps=continuation_steps,
             guess_generator=guess_maker)
