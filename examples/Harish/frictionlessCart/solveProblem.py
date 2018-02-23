from math import *
import numpy
from odesAndBcs import *
import beluga.bvpsol.BVP as BVP
from beluga.utils import ode45
import beluga.bvpsol.algorithms.MultipleShooting as MultipleShooting

x0 = 0
v0 = 0
lamx0 = 0.1
lamv0 = -0.1
tf0 = 0.1

x0 = numpy.array([x0,v0,lamx0,lamv0,tf0])
param_guess = numpy.array([])
xCont = 1
vCont = 0
flag = 1

aux = [xCont, vCont, flag]

bvp = BVP(odes,bcs,dae_func_gen=[],dae_num_states=0)
bvp.solution.aux = aux


tspan = [0, 1]

#print(bcs(x0,x0,params,aux))
[t,x] = ode45(bvp.deriv_func,tspan,x0,param_guess,bvp.solution.aux)
#[t,x] = ode45(odes,tspan,x0,param_guess,hCont, sCont, epsi, flag)

bvp.solution.x = t
bvp.solution.y = x.T
bvp.solution.parameters = param_guess

n = 100
xConts = numpy.linspace(x[-1][0],1,n)
vConts = numpy.linspace(x[-1][1],0,n)


solver = MultipleShooting(tolerance=1e-3, max_iterations=10000, max_error=10000, derivative_method='fd', cache_dir = None,verbose=True,cached=True,number_arcs=1)

for idx in range(n):
    bvp.solution.aux[0] = xConts[idx]
    bvp.solution.aux[1] = vConts[idx]
    sol = solver.solve(bvp)
    bvp.solution = sol
    print(idx)

n=100
aux[2] = 2
vConts = numpy.linspace(bvp.solution.y[1][-1],0,n)
for idx in range(n):
    bvp.solution.aux[0] = xConts[-1]
    bvp.solution.aux[1] = vConts[idx]
    sol = solver.solve(bvp)
    bvp.solution = sol
    print(idx)


#XDots = odes(tau,states,params,hCont, sCont, epsi, flag)
#print(x)

#bcVec = odes(x0,x0,param_guess,aux)
#print(bcVec)


#y0g = [bvp.solution.y[:,int(numpy.floor(i/2*x.shape[0]))] for i in range(2)]
#print(y0g)
