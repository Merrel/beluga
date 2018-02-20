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
theta0 = 0
w0 = 0
c0 = 4
lamv0 = -0.1
lamgam0 = 0
lamh0 = 0.1
lams0 = -0.1
lamtheta0 = 0.1
lamw0 = -0.1
lamc0 = -0.1

tf = 0.1

x0 = numpy.array([v0,gam0,h0,s0,theta0,w0,c0,lamv0,lamgam0,lamh0,lamtheta0,lamw0,lamc0])
param_guess = numpy.array([tf,lams0])
hCont = 100
sCont = 10
flag = 1
epsi = 1

aux = [hCont, sCont, epsi, flag]

bvp = BVP(odes,bcs,dae_func_gen=[],dae_num_states=[])
bvp.solution.aux = aux


tspan = [0, 1]

#print(bcs(x0,x0,params,aux))

[t,x] = ode45(bvp.deriv_func,tspan,x0,param_guess,bvp.solution.aux)
#[t,x] = ode45(odes,tspan,x0,param_guess,hCont, sCont, epsi, flag)

bvp.solution.x = t
bvp.solution.y = x.T
bvp.solution.parameters = param_guess

n = 100
hConts = numpy.linspace(x[-1][2],0,n)
sConts = numpy.linspace(x[-1][3],250/6378.137,n)

solver = MultipleShooting(tolerance=1e-3, max_iterations=100, max_error=100, derivative_method='fd', cache_dir = None,verbose=True,cached=True,number_arcs=1)

for idx in range(1):
    bvp.solution.aux[0] = hConts[idx]
    bvp.solution.aux[1] = sConts[idx]
    sol = solver.solve(bvp)
    bvp.solution = sol
    print(idx)

#XDots = odes(tau,states,params,hCont, sCont, epsi, flag)
#print(x)


#y0g = [bvp.solution.y[:,int(numpy.floor(i/2*x.shape[0]))] for i in range(2)]
#print(y0g)
