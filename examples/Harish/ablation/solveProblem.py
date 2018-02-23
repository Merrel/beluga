from math import *
import numpy
from odesAndBcs import *
import beluga.bvpsol.BVP as BVP
from beluga.utils import ode45
import beluga.bvpsol.algorithms.MultipleShooting as MultipleShooting

v0 = 1
gam0 = 0
h0 = 1
s0 = 0
c0 = 4
lamv0 = -0.1
lamgam0 = -0
lamh0 = 0.1
lams0 = -0.1
lamc0 = -0.1

tf = 0.1

x0 = numpy.array([v0,gam0,h0,s0,c0,lamv0,lamgam0,lamh0,lams0,lamc0,tf])

param_guess = numpy.array([])

hCont = 100.0
sCont = 10.0
flag = 1

aux = [hCont, sCont, flag]

bvp = BVP(odes,bcs,dae_func_gen=[],dae_num_states=0)
bvp.solution.aux = aux

tspan = [0, 1]
print('Starting initial guess')
[t,x] = ode45(bvp.deriv_func,tspan,x0,param_guess,bvp.solution.aux)
print('Initial guess generated')

bvp.solution.x = t
bvp.solution.y = x.T
bvp.solution.parameters = param_guess

n = 100

hConts = numpy.linspace(x[-1][2],0,n)
sConts = numpy.linspace(x[-1][3],540/6378.137,n)

solver = MultipleShooting(tolerance=1e-3, max_iterations=10000, max_error=10000, derivative_method='fd', cache_dir = None,verbose=True,cached=True,number_arcs=1)


for idx in range(n):
    bvp.solution.aux[0] = hConts[idx]
    bvp.solution.aux[1] = sConts[idx]
    sol = solver.solve(bvp)
    bvp.solution = sol
    print(idx)

n = 100
hConts = numpy.linspace(bvp.solution.y[2][-1],0,n)
bvp.solution.aux[2] = 2

for idx in range(n):
    bvp.solution.aux[0] = hConts[idx]
    bvp.solution.aux[1] = sConts[-1]
    sol = solver.solve(bvp)
    bvp.solution = sol
    print(idx)
