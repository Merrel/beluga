from math import *
import numpy
from scipy.optimize import fsolve
from math import exp as exp2

def exp(x):
    if x>200:
        x = 200
    elif x<-200:
        x = -200
    return exp2(x)

def odes(tau, states, params, aux):

    hCont = aux[0]
    sCont = aux[1]
    epsi = aux[2]
    flag = aux[3]

    v = states[0]
    gam = states[1]
    h = states[2]
    s = states[3]
    theta = states[4]
    w = states[5]
    c = states[6]

    lamv = states[7]
    lamgam = states[8]
    lamh = states[9]
    lams = states[10]
    lamtheta = states[11]
    lamw = states[12]
    lamc = states[13]

    tf = states[14]
    t = tau*tf

    mu = 398600.4418e9
    R = 6378137.0

    k = 1.74153e-4
    QStar = (8.5883e6)/3

    rho0 = 1.225
    HScale = 8500.0

    h0=40000.0
    v0=3000.0

    # Fuselage
    rBaseFuselage = 0.5
    densityFuselage = 800.0

    # Warhead
    densityWarhead = 17000.0
    rWarhead = 0.2
    warheadPosition = 2.5

    # Subsystem
    lengthSubsystem = 2.2
    rBaseSubsystem = 0.4
    densitySubsystem = 100

    #u = bisect(controlFun, -1e5, 1e5, args=(states, t, epsi), xtol=1e-3, rtol=1e-3, maxiter=10000, full_output=False, disp=False)
    u = fsolve(controlFun, 0, args=(states, t, epsi), fprime=None, full_output=0, col_deriv=0, xtol=1e-3, maxfev=0, band=None, epsfcn=None, factor=100, diag=None)

    aeroCoeff = aero_autodiff(theta,gam,v,h,w,c,u)

    fx = aeroCoeff[0]
    fz = aeroCoeff[1]
    My = aeroCoeff[2]

    fxdot_h = aeroCoeff[3]
    fzdot_h = aeroCoeff[4]
    Mydot_h = aeroCoeff[5]

    fxdot_v = aeroCoeff[6]
    fzdot_v = aeroCoeff[7]
    Mydot_v = aeroCoeff[8]

    fxdot_gam = aeroCoeff[9]
    fzdot_gam = aeroCoeff[10]
    Mydot_gam = aeroCoeff[11]

    fxdot_theta = aeroCoeff[12]
    fzdot_theta = aeroCoeff[13]
    Mydot_theta = aeroCoeff[14]

    fxdot_w = aeroCoeff[15]
    fzdot_w = aeroCoeff[16]
    Mydot_w = aeroCoeff[17]

    fxdot_c = aeroCoeff[18]
    fzdot_c = aeroCoeff[19]
    Mydot_c = aeroCoeff[20]


    vdot = ((cos(gam - theta)*fx - sin(gam - theta)*fz)/((4*densityWarhead*rWarhead**3*pi)/3 \
     - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 \
      - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + \
      (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2) - (mu*sin(gam))/(R + h*h0)**2)/v0

    gamdot = (v*v0*cos(gam))/(R + h*h0) - (cos(gam - theta)*fz + sin(gam - theta)*fx) \
    /(v*v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 \
    + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 \
    + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) - (mu*cos(gam))/(v*v0*(R + h*h0)**2)

    hdot = (v*v0*sin(gam))/h0

    sdot = (v*v0*cos(gam))/(R + h*h0)

    thetadot = w

    wdot = (cos(gam)*((cos(gam - theta)*fx - sin(gam - theta)*fz)/((4*densityWarhead*rWarhead**3*pi)/3 \
    - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 \
    - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + \
    (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2) - \
    (mu*sin(gam))/(R + h*h0)**2))/(R + h*h0) - (My + (2**(1/2)*k*v**3*v0**3*(w - \
    (v*v0*cos(gam))/(R + h*h0))*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2)\
    *((c**2*densityFuselage*rBaseFuselage**2*pi)/6 + (densityFuselage*rBaseFuselage**2*pi*\
    (c**2 + rBaseFuselage**2))/12 + (densityFuselage*rBaseFuselage**2*pi*\
    ((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3\
    *warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - \
    (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + \
    (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2)/(2*\
    ((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + \
    (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2\
    + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2) - \
    (2*c*densityFuselage*rBaseFuselage**2*pi*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 \
    - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + \
    (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2\
    *rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6))\
    /(3*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + \
    (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 \
    + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2))))/(QStar*densityFuselage) \
    - (2**(1/2)*k*rBaseFuselage**2*v**3*v0**3*pi*(w - (v*v0*cos(gam))/(R + h*h0))\
    *((6*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - \
    (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + \
    (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2\
    *rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2)\
    /((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + \
    (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*\
    rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2\
    + 3*c**2 + rBaseFuselage**2 - (8*c*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3\
    - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + \
    (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2\
    *rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6))\
    /((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 \
    + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem\
    *rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2))\
    *((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2))/(12*QStar))/\
    (((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*\
    rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 \
    - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + \
    (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2/((4*densityWarhead*rWarhead**3*pi)/3\
    - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 \
    - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + \
    (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2) - \
    (8*densityWarhead*rWarhead**5*pi)/15 + (8*densityFuselage*rWarhead**5*pi)/15 \
    - (4*densityWarhead*rWarhead**3*warheadPosition**2*pi)/3 + \
    (4*densityFuselage*rWarhead**3*warheadPosition**2*pi)/3 - \
    (c*densityFuselage*rBaseFuselage**2*pi*(c**2 + rBaseFuselage**2))/12 + \
    (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi*(lengthSubsystem**2 \
    + rBaseSubsystem**2))/12 - (densitySubsystem*lengthSubsystem*rBaseSubsystem**2\
    *pi*(lengthSubsystem**2 + rBaseSubsystem**2))/12) + (v*v0*sin(gam)*((cos(gam - theta)\
    *fz + sin(gam - theta)*fx)/(v*v0*((4*densityWarhead*rWarhead**3*pi)/3 - \
    (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 \
    - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + \
    (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) - \
    (v*v0*cos(gam))/(R + h*h0) + (mu*cos(gam))/(v*v0*(R + h*h0)**2)))/(R + h*h0) \
    - (v**2*v0**2*cos(gam)*sin(gam))/(R + h*h0)**2

    cdot = -(2**(1/2)*k*v**3*v0**3*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2))/(QStar*densityFuselage)

    lamvdot = lamw*((Mydot_v + (3*2**(1/2)*k*v**2*v0**3*(w - (v*v0*cos(gam))/(R + h*h0))*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2)*((c**2*densityFuselage*rBaseFuselage**2*pi)/6 + (densityFuselage*rBaseFuselage**2*pi*(c**2 + rBaseFuselage**2))/12 + (densityFuselage*rBaseFuselage**2*pi*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2)/(2*((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2) - (2*c*densityFuselage*rBaseFuselage**2*pi*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6))/(3*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2))))/(QStar*densityFuselage) - (2**(1/2)*k*rBaseFuselage**2*v**2*v0**3*pi*(w - (v*v0*cos(gam))/(R + h*h0))*((6*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2)/((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2 + 3*c**2 + rBaseFuselage**2 - (8*c*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6))/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2))*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2))/(4*QStar) - (2**(1/2)*k*v**3*v0**4*cos(gam)*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2)*((c**2*densityFuselage*rBaseFuselage**2*pi)/6 + (densityFuselage*rBaseFuselage**2*pi*(c**2 + rBaseFuselage**2))/12 + (densityFuselage*rBaseFuselage**2*pi*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2)/(2*((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2) - (2*c*densityFuselage*rBaseFuselage**2*pi*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6))/(3*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2))))/(QStar*densityFuselage*(R + h*h0)) + (2**(1/2)*k*rBaseFuselage**2*v**3*v0**4*pi*cos(gam)*((6*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2)/((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2 + 3*c**2 + rBaseFuselage**2 - (8*c*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6))/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2))*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2))/(12*QStar*(R + h*h0)))/(((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2) - (8*densityWarhead*rWarhead**5*pi)/15 + (8*densityFuselage*rWarhead**5*pi)/15 - (4*densityWarhead*rWarhead**3*warheadPosition**2*pi)/3 + (4*densityFuselage*rWarhead**3*warheadPosition**2*pi)/3 - (c*densityFuselage*rBaseFuselage**2*pi*(c**2 + rBaseFuselage**2))/12 + (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi*(lengthSubsystem**2 + rBaseSubsystem**2))/12 - (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi*(lengthSubsystem**2 + rBaseSubsystem**2))/12) - (cos(gam)*(cos(gam - theta)*fxdot_v - sin(gam - theta)*fzdot_v))/((R + h*h0)*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) - (v0*sin(gam)*((cos(gam - theta)*fz + sin(gam - theta)*fx)/(v*v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) - (v*v0*cos(gam))/(R + h*h0) + (mu*cos(gam))/(v*v0*(R + h*h0)**2)))/(R + h*h0) + (2*v*v0**2*cos(gam)*sin(gam))/(R + h*h0)**2 + (v*v0*sin(gam)*((v0*cos(gam))/(R + h*h0) - (cos(gam - theta)*fzdot_v + sin(gam - theta)*fxdot_v)/(v*v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) + (cos(gam - theta)*fz + sin(gam - theta)*fx)/(v**2*v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) + (mu*cos(gam))/(v**2*v0*(R + h*h0)**2)))/(R + h*h0)) - lamgam*((v0*cos(gam))/(R + h*h0) - (cos(gam - theta)*fzdot_v + sin(gam - theta)*fxdot_v)/(v*v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) + (cos(gam - theta)*fz + sin(gam - theta)*fx)/(v**2*v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) + (mu*cos(gam))/(v**2*v0*(R + h*h0)**2)) - (lams*v0*cos(gam))/(R + h*h0) - (lamv*(cos(gam - theta)*fxdot_v - sin(gam - theta)*fzdot_v))/(v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) - (lamh*v0*sin(gam))/h0 + (3*2**(1/2)*k*lamc*v**2*v0**3*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2))/(QStar*densityFuselage)

    lamgamdot = lamw*((Mydot_gam + (2**(1/2)*k*v**4*v0**4*sin(gam)*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2)*((c**2*densityFuselage*rBaseFuselage**2*pi)/6 + (densityFuselage*rBaseFuselage**2*pi*(c**2 + rBaseFuselage**2))/12 + (densityFuselage*rBaseFuselage**2*pi*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2)/(2*((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2) - (2*c*densityFuselage*rBaseFuselage**2*pi*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6))/(3*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2))))/(QStar*densityFuselage*(R + h*h0)) - (2**(1/2)*k*rBaseFuselage**2*v**4*v0**4*pi*sin(gam)*((6*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2)/((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2 + 3*c**2 + rBaseFuselage**2 - (8*c*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6))/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2))*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2))/(12*QStar*(R + h*h0)))/(((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2) - (8*densityWarhead*rWarhead**5*pi)/15 + (8*densityFuselage*rWarhead**5*pi)/15 - (4*densityWarhead*rWarhead**3*warheadPosition**2*pi)/3 + (4*densityFuselage*rWarhead**3*warheadPosition**2*pi)/3 - (c*densityFuselage*rBaseFuselage**2*pi*(c**2 + rBaseFuselage**2))/12 + (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi*(lengthSubsystem**2 + rBaseSubsystem**2))/12 - (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi*(lengthSubsystem**2 + rBaseSubsystem**2))/12) + (sin(gam)*((cos(gam - theta)*fx - sin(gam - theta)*fz)/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2) - (mu*sin(gam))/(R + h*h0)**2))/(R + h*h0) + (cos(gam)*((sin(gam - theta)*fzdot_gam - cos(gam - theta)*fxdot_gam + cos(gam - theta)*fz + sin(gam - theta)*fx)/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2) + (mu*cos(gam))/(R + h*h0)**2))/(R + h*h0) + (v**2*v0**2*cos(gam)**2)/(R + h*h0)**2 - (v**2*v0**2*sin(gam)**2)/(R + h*h0)**2 - (v*v0*cos(gam)*((cos(gam - theta)*fz + sin(gam - theta)*fx)/(v*v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) - (v*v0*cos(gam))/(R + h*h0) + (mu*cos(gam))/(v*v0*(R + h*h0)**2)))/(R + h*h0) - (v*v0*sin(gam)*((v*v0*sin(gam))/(R + h*h0) + (cos(gam - theta)*fzdot_gam + sin(gam - theta)*fxdot_gam + cos(gam - theta)*fx - sin(gam - theta)*fz)/(v*v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) - (mu*sin(gam))/(v*v0*(R + h*h0)**2)))/(R + h*h0)) + lamgam*((v*v0*sin(gam))/(R + h*h0) + (cos(gam - theta)*fzdot_gam + sin(gam - theta)*fxdot_gam + cos(gam - theta)*fx - sin(gam - theta)*fz)/(v*v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) - (mu*sin(gam))/(v*v0*(R + h*h0)**2)) + (lamv*((sin(gam - theta)*fzdot_gam - cos(gam - theta)*fxdot_gam + cos(gam - theta)*fz + sin(gam - theta)*fx)/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2) + (mu*cos(gam))/(R + h*h0)**2))/v0 - (lamh*v*v0*cos(gam))/h0 + (lams*v*v0*sin(gam))/(R + h*h0)

    lamhdot = lamw*((Mydot_h + (2**(1/2)*h0*k*v**4*v0**4*cos(gam)*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2)*((c**2*densityFuselage*rBaseFuselage**2*pi)/6 + (densityFuselage*rBaseFuselage**2*pi*(c**2 + rBaseFuselage**2))/12 + (densityFuselage*rBaseFuselage**2*pi*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2)/(2*((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2) - (2*c*densityFuselage*rBaseFuselage**2*pi*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6))/(3*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2))))/(QStar*densityFuselage*(R + h*h0)**2) - (2**(1/2)*h0*k*rBaseFuselage**2*v**4*v0**4*pi*cos(gam)*((6*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2)/((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2 + 3*c**2 + rBaseFuselage**2 - (8*c*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6))/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2))*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2))/(12*QStar*(R + h*h0)**2) + (2**(1/2)*c*h0*k*rho0*v**3*v0**3*pi*exp(-(h*h0)/HScale)*(w - (v*v0*cos(gam))/(R + h*h0))*((6*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2)/((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2 + 3*c**2 + rBaseFuselage**2 - (8*c*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6))/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)))/(24*HScale*QStar*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2)) - (2**(1/2)*c*h0*k*rho0*v**3*v0**3*exp(-(h*h0)/HScale)*(w - (v*v0*cos(gam))/(R + h*h0))*((c**2*densityFuselage*rBaseFuselage**2*pi)/6 + (densityFuselage*rBaseFuselage**2*pi*(c**2 + rBaseFuselage**2))/12 + (densityFuselage*rBaseFuselage**2*pi*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2)/(2*((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2) - (2*c*densityFuselage*rBaseFuselage**2*pi*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6))/(3*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2))))/(2*HScale*QStar*densityFuselage*rBaseFuselage**2*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2)))/(((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2) - (8*densityWarhead*rWarhead**5*pi)/15 + (8*densityFuselage*rWarhead**5*pi)/15 - (4*densityWarhead*rWarhead**3*warheadPosition**2*pi)/3 + (4*densityFuselage*rWarhead**3*warheadPosition**2*pi)/3 - (c*densityFuselage*rBaseFuselage**2*pi*(c**2 + rBaseFuselage**2))/12 + (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi*(lengthSubsystem**2 + rBaseSubsystem**2))/12 - (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi*(lengthSubsystem**2 + rBaseSubsystem**2))/12) - (cos(gam)*((cos(gam - theta)*fxdot_h - sin(gam - theta)*fzdot_h)/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2) + (2*h0*mu*sin(gam))/(R + h*h0)**3))/(R + h*h0) + (h0*cos(gam)*((cos(gam - theta)*fx - sin(gam - theta)*fz)/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2) - (mu*sin(gam))/(R + h*h0)**2))/(R + h*h0)**2 - (v*v0*sin(gam)*((cos(gam - theta)*fzdot_h + sin(gam - theta)*fxdot_h)/(v*v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) + (h0*v*v0*cos(gam))/(R + h*h0)**2 - (2*h0*mu*cos(gam))/(v*v0*(R + h*h0)**3)))/(R + h*h0) - (2*h0*v**2*v0**2*cos(gam)*sin(gam))/(R + h*h0)**3 + (h0*v*v0*sin(gam)*((cos(gam - theta)*fz + sin(gam - theta)*fx)/(v*v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) - (v*v0*cos(gam))/(R + h*h0) + (mu*cos(gam))/(v*v0*(R + h*h0)**2)))/(R + h*h0)**2) + lamgam*((cos(gam - theta)*fzdot_h + sin(gam - theta)*fxdot_h)/(v*v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) + (h0*v*v0*cos(gam))/(R + h*h0)**2 - (2*h0*mu*cos(gam))/(v*v0*(R + h*h0)**3)) - (lamv*((cos(gam - theta)*fxdot_h - sin(gam - theta)*fzdot_h)/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2) + (2*h0*mu*sin(gam))/(R + h*h0)**3))/v0 + (h0*lams*v*v0*cos(gam))/(R + h*h0)**2 - (2**(1/2)*c*h0*k*lamc*rho0*v**3*v0**3*exp(-(h*h0)/HScale))/(2*HScale*QStar*densityFuselage*rBaseFuselage**2*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2))

    lamsdot = 0

    lamthetadot = (lamgam*(cos(gam - theta)*fzdot_theta + sin(gam - theta)*fxdot_theta - cos(gam - theta)*fx + sin(gam - theta)*fz))/(v*v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) - (lamv*(cos(gam - theta)*fxdot_theta - sin(gam - theta)*fzdot_theta + cos(gam - theta)*fz + sin(gam - theta)*fx))/(v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) - lamw*((cos(gam)*(cos(gam - theta)*fxdot_theta - sin(gam - theta)*fzdot_theta + cos(gam - theta)*fz + sin(gam - theta)*fx))/((R + h*h0)*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) - Mydot_theta/(((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2) - (8*densityWarhead*rWarhead**5*pi)/15 + (8*densityFuselage*rWarhead**5*pi)/15 - (4*densityWarhead*rWarhead**3*warheadPosition**2*pi)/3 + (4*densityFuselage*rWarhead**3*warheadPosition**2*pi)/3 - (c*densityFuselage*rBaseFuselage**2*pi*(c**2 + rBaseFuselage**2))/12 + (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi*(lengthSubsystem**2 + rBaseSubsystem**2))/12 - (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi*(lengthSubsystem**2 + rBaseSubsystem**2))/12) + (sin(gam)*(cos(gam - theta)*fzdot_theta + sin(gam - theta)*fxdot_theta - cos(gam - theta)*fx + sin(gam - theta)*fz))/((R + h*h0)*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)))

    lamwdot = (lamgam*(cos(gam - theta)*fzdot_w + sin(gam - theta)*fxdot_w))/(v*v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) - lamw*((cos(gam)*(cos(gam - theta)*fxdot_w - sin(gam - theta)*fzdot_w))/((R + h*h0)*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) - (Mydot_w + (2**(1/2)*k*v**3*v0**3*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2)*((c**2*densityFuselage*rBaseFuselage**2*pi)/6 + (densityFuselage*rBaseFuselage**2*pi*(c**2 + rBaseFuselage**2))/12 + (densityFuselage*rBaseFuselage**2*pi*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2)/(2*((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2) - (2*c*densityFuselage*rBaseFuselage**2*pi*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6))/(3*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2))))/(QStar*densityFuselage) - (2**(1/2)*k*rBaseFuselage**2*v**3*v0**3*pi*((6*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2)/((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2 + 3*c**2 + rBaseFuselage**2 - (8*c*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6))/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2))*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2))/(12*QStar))/(((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2) - (8*densityWarhead*rWarhead**5*pi)/15 + (8*densityFuselage*rWarhead**5*pi)/15 - (4*densityWarhead*rWarhead**3*warheadPosition**2*pi)/3 + (4*densityFuselage*rWarhead**3*warheadPosition**2*pi)/3 - (c*densityFuselage*rBaseFuselage**2*pi*(c**2 + rBaseFuselage**2))/12 + (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi*(lengthSubsystem**2 + rBaseSubsystem**2))/12 - (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi*(lengthSubsystem**2 + rBaseSubsystem**2))/12) + (sin(gam)*(cos(gam - theta)*fzdot_w + sin(gam - theta)*fxdot_w))/((R + h*h0)*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2))) - (lamv*(cos(gam - theta)*fxdot_w - sin(gam - theta)*fzdot_w))/(v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) - lamtheta

    lamcdot = lamgam*((cos(gam - theta)*fzdot_c + sin(gam - theta)*fxdot_c)/(v*v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) - (densityFuselage*rBaseFuselage**2*pi*(cos(gam - theta)*fz + sin(gam - theta)*fx))/(2*v*v0*((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2)) + lamw*((Mydot_c - (2**(1/2)*k*v**3*v0**3*(w - (v*v0*cos(gam))/(R + h*h0))*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2)*((2*densityFuselage*rBaseFuselage**2*pi*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6))/(3*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) - (c*densityFuselage*rBaseFuselage**2*pi)/2 + (densityFuselage**2*rBaseFuselage**4*pi**2*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2)/(2*((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**3) + (2*c**2*densityFuselage**2*rBaseFuselage**4*pi**2)/(9*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) - (2*c*densityFuselage**2*rBaseFuselage**4*pi**2*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6))/(3*((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2)))/(QStar*densityFuselage) + (2**(1/2)*k*rBaseFuselage**2*v**3*v0**3*pi*(w - (v*v0*cos(gam))/(R + h*h0))*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2)*((8*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6))/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2) - 6*c + (8*c**2*densityFuselage*rBaseFuselage**2*pi)/(3*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) + (6*densityFuselage*rBaseFuselage**2*pi*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2)/((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**3 - (8*c*densityFuselage*rBaseFuselage**2*pi*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6))/((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2))/(12*QStar) - (2**(1/2)*k*rho0*v**3*v0**3*pi*exp(-(h*h0)/HScale)*(w - (v*v0*cos(gam))/(R + h*h0))*((6*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2)/((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2 + 3*c**2 + rBaseFuselage**2 - (8*c*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6))/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)))/(24*QStar*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2)) + (2**(1/2)*k*rho0*v**3*v0**3*exp(-(h*h0)/HScale)*(w - (v*v0*cos(gam))/(R + h*h0))*((c**2*densityFuselage*rBaseFuselage**2*pi)/6 + (densityFuselage*rBaseFuselage**2*pi*(c**2 + rBaseFuselage**2))/12 + (densityFuselage*rBaseFuselage**2*pi*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2)/(2*((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2) - (2*c*densityFuselage*rBaseFuselage**2*pi*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6))/(3*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2))))/(2*QStar*densityFuselage*rBaseFuselage**2*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2)))/(((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2) - (8*densityWarhead*rWarhead**5*pi)/15 + (8*densityFuselage*rWarhead**5*pi)/15 - (4*densityWarhead*rWarhead**3*warheadPosition**2*pi)/3 + (4*densityFuselage*rWarhead**3*warheadPosition**2*pi)/3 - (c*densityFuselage*rBaseFuselage**2*pi*(c**2 + rBaseFuselage**2))/12 + (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi*(lengthSubsystem**2 + rBaseSubsystem**2))/12 - (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi*(lengthSubsystem**2 + rBaseSubsystem**2))/12) + ((My + (2**(1/2)*k*v**3*v0**3*(w - (v*v0*cos(gam))/(R + h*h0))*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2)*((c**2*densityFuselage*rBaseFuselage**2*pi)/6 + (densityFuselage*rBaseFuselage**2*pi*(c**2 + rBaseFuselage**2))/12 + (densityFuselage*rBaseFuselage**2*pi*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2)/(2*((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2) - (2*c*densityFuselage*rBaseFuselage**2*pi*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6))/(3*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2))))/(QStar*densityFuselage) - (2**(1/2)*k*rBaseFuselage**2*v**3*v0**3*pi*(w - (v*v0*cos(gam))/(R + h*h0))*((6*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2)/((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2 + 3*c**2 + rBaseFuselage**2 - (8*c*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6))/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2))*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2))/(12*QStar))*((c**2*densityFuselage*rBaseFuselage**2*pi)/6 + (densityFuselage*rBaseFuselage**2*pi*(c**2 + rBaseFuselage**2))/12 + (densityFuselage*rBaseFuselage**2*pi*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2)/(2*((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2) - (2*c*densityFuselage*rBaseFuselage**2*pi*((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6))/(3*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2))))/(((4*densityWarhead*pi*rWarhead**3*warheadPosition)/3 - (4*densityFuselage*pi*rWarhead**3*warheadPosition)/3 + (c**2*densityFuselage*pi*rBaseFuselage**2)/6 - (densityFuselage*lengthSubsystem**2*pi*rBaseSubsystem**2)/6 + (densitySubsystem*lengthSubsystem**2*pi*rBaseSubsystem**2)/6)**2/((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2) - (8*densityWarhead*pi*rWarhead**5)/15 + (8*densityFuselage*pi*rWarhead**5)/15 - (4*densityWarhead*pi*rWarhead**3*warheadPosition**2)/3 + (4*densityFuselage*pi*rWarhead**3*warheadPosition**2)/3 - (c*densityFuselage*pi*rBaseFuselage**2*(c**2 + rBaseFuselage**2))/12 + (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2*(lengthSubsystem**2 + rBaseSubsystem**2))/12 - (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2*(lengthSubsystem**2 + rBaseSubsystem**2))/12)**2 - (cos(gam)*((cos(gam - theta)*fxdot_c - sin(gam - theta)*fzdot_c)/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2) - (densityFuselage*rBaseFuselage**2*pi*(cos(gam - theta)*fx - sin(gam - theta)*fz))/(2*((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2)))/(R + h*h0) - (v*v0*sin(gam)*((cos(gam - theta)*fzdot_c + sin(gam - theta)*fxdot_c)/(v*v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) - (densityFuselage*rBaseFuselage**2*pi*(cos(gam - theta)*fz + sin(gam - theta)*fx))/(2*v*v0*((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2)))/(R + h*h0)) - (lamv*((cos(gam - theta)*fxdot_c - sin(gam - theta)*fzdot_c)/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2) - (densityFuselage*rBaseFuselage**2*pi*(cos(gam - theta)*fx - sin(gam - theta)*fz))/(2*((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2)))/v0 + (2**(1/2)*k*lamc*rho0*v**3*v0**3*exp(-(h*h0)/HScale))/(2*QStar*densityFuselage*rBaseFuselage**2*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2))

    tfdot = 0

    XDot = tf*numpy.array([vdot,gamdot,hdot,sdot,thetadot,wdot,cdot,lamvdot,lamgamdot,lamhdot,lamsdot,lamthetadot,lamwdot,lamcdot,tfdot])
    return XDot

def bcs(Ya, Yb, params, aux):
    hCont = aux[0]
    sCont = aux[1]
    epsi = aux[2]
    flag = aux[3]

    v0 = Ya[0]
    gam0 = Ya[1]
    h0 = Ya[2]
    s0 = Ya[3]
    theta0 = Ya[4]
    w0 = Ya[5]
    c0 = Ya[6]
    lamv0 = Ya[7]
    lamgam0 = Ya[8]
    lamh0 = Ya[9]
    lams0 = Ya[10]
    lamtheta0 = Ya[11]
    lamw0 = Ya[12]
    lamc0 = Ya[13]
    tf0 = Ya[14]

    v = Yb[0]
    gam = Yb[1]
    h = Yb[2]
    s = Yb[3]
    theta = Yb[4]
    w = Yb[5]
    c = Yb[6]
    lamv = Yb[7]
    lamgam = Yb[8]
    lamh = Yb[9]
    lams = Yb[10]
    lamtheta = Yb[11]
    lamw = Yb[12]
    lamc = Yb[13]
    tf = Yb[14]

    Xdot = odes(1, Yb, params, aux)

    #u = bisect(controlFun, -1e5, 1e5, args=(Yb, tf, epsi), xtol=1e-3, rtol=1e-3, maxiter=10000, full_output=False, disp=False)
    u = fsolve(controlFun, 0, args=(Yb, tf, epsi), fprime=None, full_output=0, col_deriv=0, xtol=1e-3, maxfev=0, band=None, epsfcn=None, factor=100, diag=None)


    weight = 0;

    coState = numpy.array([lamv, lamgam, lamh, lams, lamtheta, lamw, lamc])

    Ham = (epsi*(u**2)) + (weight*1) + (numpy.dot(coState,Xdot[:7])/tf)

    if flag == 1:
        zeroVec = numpy.array([v0 - 1,\
        gam0 - 0,\
        h0 - 1,\
        s0 - 0,\
        theta0 - 0,\
        w0 - 0,\
        c0 - 4,\
        lamv + v,\
        lamgam - 0,\
        lamh - 0,\
        s - sCont,\
        lamtheta - 0,\
        lamw - 0,\
        lamc - 0,\
        Ham - 0])
    else:
        zeroVec = numpy.array([v0 - 1,\
        gam0 - 0,\
        h0 - 1,\
        s0 - 0,\
        theta0 - 0,\
        w0 - 0,\
        c0 - 4,\
        lamv + v,\
        lamgam - 0,\
        h - hCont,\
        s - sCont,\
        lamtheta - 0,\
        lamw - 0,\
        lamc - 0,\
        Ham - 0])

    return zeroVec

def controlFun(u,states, t, epsi):
    mu = 398600.4418e9
    R = 6378137.0

    k = 1.74153e-4
    QStar = (8.5883e6)/3

    rho0 = 1.225
    HScale = 8500

    h0=40000.0
    v0=3000.0

    # Fuselage
    rBaseFuselage = 0.5
    densityFuselage = 800

    # Warhead
    densityWarhead = 17000.0
    rWarhead = 0.2
    warheadPosition = 2.5

    # Subsystem
    lengthSubsystem = 2.2
    rBaseSubsystem = 0.4
    densitySubsystem = 100

    v = states[0]
    gam = states[1]
    h = states[2]
    s = states[3]
    theta = states[4]
    w = states[5]
    c = states[6]

    lamv = states[7]
    lamgam = states[8]
    lamh = states[9]
    lams = states[10]
    lamtheta = states[11]
    lamw = states[12]
    lamc = states[13]

    tf = states[14]

    aeroCoeff = aero_autodiff_u(theta,gam,v,h,w,c,u)
    fxdot_u = aeroCoeff[0]
    fzdot_u = aeroCoeff[1]
    Mydot_u = aeroCoeff[2]

    dHdu = 2*epsi*u + lamw*((cos(gam)*(cos(gam - theta)*fxdot_u - sin(gam - theta)*fzdot_u))/((R + h*h0)*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) - Mydot_u/(((4*densityWarhead*rWarhead**3*warheadPosition*pi)/3 - (4*densityFuselage*rWarhead**3*warheadPosition*pi)/3 + (c**2*densityFuselage*rBaseFuselage**2*pi)/6 - (densityFuselage*lengthSubsystem**2*rBaseSubsystem**2*pi)/6 + (densitySubsystem*lengthSubsystem**2*rBaseSubsystem**2*pi)/6)**2/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2) - (8*densityWarhead*rWarhead**5*pi)/15 + (8*densityFuselage*rWarhead**5*pi)/15 - (4*densityWarhead*rWarhead**3*warheadPosition**2*pi)/3 + (4*densityFuselage*rWarhead**3*warheadPosition**2*pi)/3 - (c*densityFuselage*rBaseFuselage**2*pi*(c**2 + rBaseFuselage**2))/12 + (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi*(lengthSubsystem**2 + rBaseSubsystem**2))/12 - (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi*(lengthSubsystem**2 + rBaseSubsystem**2))/12) + (sin(gam)*(cos(gam - theta)*fzdot_u + sin(gam - theta)*fxdot_u))/((R + h*h0)*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2))) + (lamv*(cos(gam - theta)*fxdot_u - sin(gam - theta)*fzdot_u))/(v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) - (lamgam*(cos(gam - theta)*fzdot_u + sin(gam - theta)*fxdot_u))/(v*v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2))

    return dHdu

def aero_autodiff(theta,gam,V,Y,w,c,u_sat):
    '''
    This function models the hypersonic aerodynamics for an elliptic paraboloid
    vehicle.

    The aerodynamic force and moment components are computed using
    panel methods and Modified Newtonian flow theory. Unsteady flow is modeled
    using piston theory.

    The derivatives with respect to the inputs are computed using automatic
    differentiation. The necessary force and moment components and their
    derivatives are returned as a numpy array.

    The derivatives with respect to c are required when the vehicle ablates.

    Inputs:
    theta = downrange angle (longitude), rad
    gam = flight-path-angle, rad
    V = scaled velocity, nd
    Y = scaled altitude, nd
    w = pitch rate, rad/s
    c = length of fuselage, m
    u_sat = fictitious control used for bounding original control,
            elevon deflection command (del_c)
    '''

    # Compute velocity and altitude from scale factors
    v0 = 3000.0
    y0 = 40000.0
    v = V*v0
    vdot_V = v0

    # Compute control command from fictitious control
    ub_del_c=90*pi/180
    lb_del_c=-ub_del_c
    slp_u=0.3

    del_c=(ub_del_c-((ub_del_c-lb_del_c)/(1+exp(slp_u*u_sat))));

    # Compute elevon deflection angle after the influence of the washout filter
    k=1.0 # Washout folter gain
    KK=1.0

    delta=del_c+(k*w)
    deldot_w=k;

    # Compute angle-of-attack
    alpha=theta-gam
    alphadot_theta=1
    alphadot_gam=-1

    # Compute vehicle properties
    a_b = 0.5 # semi-major axis of fuselage base
    b_b = 0.5 # semi-minor axis of fuselage base

    # Fuselage properties
    rBaseFuselage = 0.5
    densityFuselage = 800.0

    # Warhead properties
    densityWarhead = 17000.0
    rWarhead = 0.2
    warheadPosition = 2.5

    # Subsystem properties
    lengthSubsystem = 2.2
    rBaseSubsystem = 0.4
    densitySubsystem = 100

    # Mass of components
    massFuselage = (densityFuselage*pi*(rBaseFuselage**2)*c/2) - \
    (densityFuselage*(4/3)*pi*(rWarhead**3)) - (densityFuselage*pi*(rBaseSubsystem**2)*lengthSubsystem/2)

    massFuselagedot_c = (densityFuselage*pi*(rBaseFuselage**2)/2)
    massWarhead = densityWarhead*(4/3)*pi*(rWarhead**3)
    massSubsystem = (densitySubsystem*pi*(rBaseSubsystem**2)*lengthSubsystem/2)
    mass = massFuselage + massWarhead + massSubsystem
    massdot_c = massFuselagedot_c

    # Center of mass of the components
    cgFuselage = (((densityFuselage*pi*(rBaseFuselage**2)*c/2)*c/3) \
    - (densityFuselage*(4/3)*pi*(rWarhead**3)*warheadPosition) \
    - ((densityFuselage*pi*(rBaseSubsystem**2)*lengthSubsystem/2)*lengthSubsystem/3))/massFuselage

    cgFuselagedot_c = (((densityFuselage*pi*(rBaseFuselage**2))*c/3)/massFuselage) \
    - ((((densityFuselage*pi*(rBaseFuselage**2)*c/2)*c/3) \
    - (densityFuselage*(4/3)*pi*(rWarhead**3)*warheadPosition) \
    - ((densityFuselage*pi*(rBaseSubsystem**2)*lengthSubsystem/2)*lengthSubsystem/3))*massFuselagedot_c/(massFuselage**2))

    cgWarhead = warheadPosition
    cgSubsystem = lengthSubsystem/3
    cg = ((massFuselage*cgFuselage) + (massWarhead*cgWarhead) + (massSubsystem*cgSubsystem))/mass
    cgdot_c = (massFuselagedot_c*cgFuselage/mass) + (massFuselage*cgFuselagedot_c/mass) \
    - (((massFuselage*cgFuselage) + (massWarhead*cgWarhead) + (massSubsystem*cgSubsystem))*massdot_c/(mass**2))

    CG = numpy.array([cg,0,0])
    CGdot_c = numpy.array([cgdot_c,0,0])

    # Freestream comditions (Replace with Earth Standard Atmosphere model in the future)

    W=numpy.array([-cos(alpha),0,-sin(alpha)]) # Relative wind
    Wdot_theta=numpy.array([sin(alpha),0,-cos(alpha)])*alphadot_theta
    Wdot_gam=numpy.array([sin(alpha),0,-cos(alpha)])*alphadot_gam

    Rgas=287
    rho0 = 1.225
    HScale = 8500
    rho=rho0*exp(-(Y*y0/HScale))
    rhodot_Y=-((rho0*y0/HScale)*exp(-(Y*y0/HScale)))

    P=101325*exp(-(Y*y0/HScale))
    Pdot_Y=-((101325*y0/HScale)*exp(-(Y*y0/HScale)))

    T=P/(rho*Rgas)
    Tdot_Y=(Pdot_Y/(rho*Rgas))-(P*rhodot_Y/((rho**2)*Rgas))

    gamma=1.4

    Mach=v/sqrt(gamma*Rgas*T)
    Machdot_V=vdot_V/sqrt(gamma*Rgas*T)
    Machdot_Y=-((1/2)*v*((gamma*Rgas*T)**(-3/2))*gamma*Rgas*Tdot_Y)

    epsi=(((gamma-1)*(Mach**2))+2)/((gamma+1)*(Mach**2))
    epsidot_Y=((2*(gamma-1)/((gamma+1)*Mach))-(2*(((gamma-1)*(Mach**2))+2)/((gamma+1)*(Mach**3))))*Machdot_Y
    epsidot_V=((2*(gamma-1)/((gamma+1)*Mach))-(2*(((gamma-1)*(Mach**2))+2)/((gamma+1)*(Mach**3))))*Machdot_V

    # Downstream conditions behind the normal portion of the shock
    P1=P*(((2*gamma*(Mach**2))-(gamma-1))/(gamma+1))
    P1dot_Y=(Pdot_Y*(((2*gamma*(Mach**2))-(gamma-1))/(gamma+1)))+(P*4*gamma*Mach*Machdot_Y/(gamma+1))
    P1dot_V=(P*4*gamma*Mach*Machdot_V/(gamma+1))

    rho1=rho/epsi
    rho1dot_Y=(rhodot_Y/epsi)-(rho*epsidot_Y/(epsi**2))
    rho1dot_V=-(rho*epsidot_V/(epsi**2))

    T1=P1/(rho1*Rgas)
    T1dot_Y=(P1dot_Y/(rho1*Rgas))-(P1*rho1dot_Y/((rho1**2)*Rgas))
    T1dot_V=(P1dot_V/(rho1*Rgas))-(P1*rho1dot_V/((rho1**2)*Rgas))

    # Initialize force vector and its derivatives
    F = numpy.array([0.0,0.0,0.0])
    Fdot_theta = numpy.array([0.0,0.0,0.0])
    Fdot_gam = numpy.array([0.0,0.0,0.0])
    Fdot_V = numpy.array([0.0,0.0,0.0])
    Fdot_Y = numpy.array([0.0,0.0,0.0])
    Fdot_w = numpy.array([0.0,0.0,0.0])
    Fdot_c = numpy.array([0.0,0.0,0.0])

    # Initialize moment vector and its derivatives
    M = numpy.array([0.0,0.0,0.0])
    Mdot_theta = numpy.array([0.0,0.0,0.0])
    Mdot_gam = numpy.array([0.0,0.0,0.0])
    Mdot_V = numpy.array([0.0,0.0,0.0])
    Mdot_Y = numpy.array([0.0,0.0,0.0])
    Mdot_w = numpy.array([0.0,0.0,0.0])
    Mdot_c = numpy.array([0.0,0.0,0.0])

    n=11 # Discretization parameter for panel methods

    # Fuselage
    x = numpy.linspace(0,c,n)
    xdot_c = numpy.linspace(0,1,n)
    thetas = numpy.linspace(pi/2,3*pi/2,n)

    for ii in range(n-1):
        for jj in range(n-1):
            if ii != n-2:
                # Compute the vertices of the quadrilateral panels
                X = x[ii]
                Xdot_c = xdot_c[ii]

                Xp=x[ii+1]
                Xpdot_c = xdot_c[ii+1]

                thetass = thetas[jj]
                thetap = thetas[jj+1]

                p1 = numpy.array([X,a_b*sqrt(1-(X/c))*cos(thetass),b_b*sqrt(1-(X/c))*sin(thetass)])
                p1dot_c = numpy.array([Xdot_c, a_b*(-(Xdot_c/c)+(X/(c**2)))*cos(thetass)/(2*sqrt(1-(X/c))), b_b*(-(Xdot_c/c)+(X/(c**2)))*sin(thetass)/(2*sqrt(1-(X/c)))])

                p2 = numpy.array([Xp, a_b*sqrt(1-(Xp/c))*cos(thetass), b_b*sqrt(1-(Xp/c))*sin(thetass)])
                p2dot_c = numpy.array([Xpdot_c, a_b*(-(Xpdot_c/c)+(Xp/(c**2)))*cos(thetass)/(2*sqrt(1-(Xp/c))), b_b*(-(Xpdot_c/c)+(Xp/(c**2)))*sin(thetass)/(2*sqrt(1-(Xp/c)))])

                p3 = numpy.array([Xp, a_b*sqrt(1-(Xp/c))*cos(thetap), b_b*sqrt(1-(Xp/c))*sin(thetap)])
                p3dot_c = numpy.array([Xpdot_c, a_b*(-(Xpdot_c/c)+(Xp/(c**2)))*cos(thetap)/(2*sqrt(1-(Xp/c))), b_b*(-(Xpdot_c/c)+(Xp/(c**2)))*sin(thetap)/(2*sqrt(1-(Xp/c)))])

                p4 = numpy.array([X, a_b*sqrt(1-(X/c))*cos(thetap), b_b*sqrt(1-(X/c))*sin(thetap)])
                p4dot_c = numpy.array([Xdot_c, a_b*(-(Xdot_c/c)+(X/(c**2)))*cos(thetap)/(2*sqrt(1-(X/c))), b_b*(-(Xdot_c/c)+(X/(c**2)))*sin(thetap)/(2*sqrt(1-(X/c)))])

                pc=(p1+p2+p3+p4)/4 # Centroid
                pcdot_c = (p1dot_c + p2dot_c + p3dot_c + p4dot_c)/4

                N = numpy.cross(p2-p1,p3-p1) # Normal vector
                Ndot_c = numpy.cross(p2dot_c-p1dot_c,p3-p1) + numpy.cross(p2-p1,p3dot_c-p1dot_c)

                normN = numpy.linalg.norm(N)
                normNdot_c = ((N[0]*Ndot_c[0]) + (N[1]*Ndot_c[1]) + (N[2]*Ndot_c[2]))/normN

                n_cap=N/numpy.linalg.norm(N) # Unit normal vector
                n_capdot_c = (Ndot_c/normN) - (N*normNdot_c/(normN**2))

                # Compute area of panel
                cross1 = numpy.cross(p2-p1,p1-p3)
                cross1dot_c = numpy.cross(p2dot_c-p1dot_c,p1-p3) + numpy.cross(p2-p1,p1dot_c-p3dot_c)
                normcross1 = numpy.linalg.norm(cross1)
                normcross1dot_c = ((cross1[0]*cross1dot_c[0]) + (cross1[1]*cross1dot_c[1]) + (cross1[2]*cross1dot_c[2]))/normcross1

                cross2 = numpy.cross(p3-p1,p1-p4)
                cross2dot_c = numpy.cross(p3dot_c-p1dot_c,p1-p4) + numpy.cross(p3-p1,p1dot_c-p4dot_c)
                normcross2 = numpy.linalg.norm(cross2)
                normcross2dot_c = ((cross2[0]*cross2dot_c[0]) + (cross2[1]*cross2dot_c[1]) + (cross2[2]*cross2dot_c[2]))/normcross2

                dA = 0.5*(normcross1 + normcross2)
                dAdot_c = 0.5*(normcross1dot_c + normcross2dot_c)
            else:
                # Compute the vertices of triangular panels
                X = x[ii];
                Xdot_c = xdot_c[ii]

                Xp = x[ii+1]
                Xpdot_c = xdot_c[ii+1]

                thetass = thetas[jj]
                thetap = thetas[jj+1]

                p1 = numpy.array([X, a_b*sqrt(1-(X/c))*cos(thetass), b_b*sqrt(1-(X/c))*sin(thetass)])
                p1dot_c = numpy.array([Xdot_c, a_b*(-(Xdot_c/c)+(X/(c**2)))*cos(thetass)/(2*sqrt(1-(X/c))), b_b*(-(Xdot_c/c)+(X/(c**2)))*sin(thetass)/(2*sqrt(1-(X/c)))])

                p2 = numpy.array([Xp, a_b*sqrt(1-(Xp/c))*cos(thetass), b_b*sqrt(1-(Xp/c))*sin(thetass)])

                if ii == n-2:
                    p2dot_c = numpy.array([Xpdot_c, 0, 0])
                else:
                    p2dot_c = numpy.array([Xpdot_c, a_b*(-(Xpdot_c/c)+(Xp/(c**2)))*cos(thetass)/(2*sqrt(1-(Xp/c))), b_b*(-(Xpdot_c/c)+(Xp/(c**2)))*sin(thetass)/(2*sqrt(1-(Xp/c)))])

                p3 = numpy.array([X, a_b*sqrt(1-(X/c))*cos(thetap), b_b*sqrt(1-(X/c))*sin(thetap)])
                p3dot_c = numpy.array([Xdot_c, a_b*(-(Xdot_c/c)+(X/(c**2)))*cos(thetap)/(2*sqrt(1-(X/c))), b_b*(-(Xdot_c/c)+(X/(c**2)))*sin(thetap)/(2*sqrt(1-(X/c)))])

                pc = (p1+p2+p3)/3 # Centroid
                pcdot_c = (p1dot_c+p2dot_c+p3dot_c)/3

                N = numpy.cross(p2-p1,p3-p1) # Normal vector
                Ndot_c = numpy.cross(p2dot_c-p1dot_c,p3-p1) + numpy.cross(p2-p1,p3dot_c-p1dot_c)

                normN = numpy.linalg.norm(N)
                normNdot_c = ((N[0]*Ndot_c[0]) + (N[1]*Ndot_c[1]) + (N[2]*Ndot_c[2]))/normN

                n_cap = N/normN # Unit normal vector
                n_capdot_c = (Ndot_c/normN) - (N*normNdot_c/(normN**2))

                # Compute area of panel
                cross1 = numpy.cross(p2-p1,p1-p3)
                cross1dot_c = numpy.cross(p2dot_c-p1dot_c,p1-p3) + numpy.cross(p2-p1,p1dot_c-p3dot_c)

                normcross1 = numpy.linalg.norm(cross1)
                normcross1dot_c = ((cross1[0]*cross1dot_c[0]) + (cross1[1]*cross1dot_c[1]) + (cross1[2]*cross1dot_c[2]))/normcross1

                dA = 0.5*normcross1
                dAdot_c = 0.5*normcross1dot_c

            # Cosine of angle between inward-pointing normal vector and relative wind
            costheta = numpy.dot(W,n_cap)
            costhetadot_theta = numpy.dot(Wdot_theta,n_cap)
            costhetadot_gam = numpy.dot(Wdot_gam,n_cap)
            costhetadot_c = numpy.dot(W,n_capdot_c)

            if costheta<0:
                # Edge conditions if panel is shadowed
                Pe=P
                Pedot_Y=Pdot_Y
                Pedot_V=0
                Pedot_gam=0
                Pedot_theta=0
                Pedot_c = 0

                Te=T
                Tedot_Y=Tdot_Y
                Tedot_V=0
                Tedot_gam=0
                Tedot_theta=0
                Tedot_c = 0

                rhoe=Pe/(Rgas*Te)
                rhoedot_Y=(Pedot_Y/(Rgas*Te))-(Pe*Tedot_Y/(Rgas*(Te**2)))
                rhoedot_V=(Pedot_V/(Rgas*Te))-(Pe*Tedot_V/(Rgas*(Te**2)))
                rhoedot_gam=(Pedot_gam/(Rgas*Te))-(Pe*Tedot_gam/(Rgas*(Te**2)))
                rhoedot_theta=(Pedot_theta/(Rgas*Te))-(Pe*Tedot_theta/(Rgas*(Te**2)))
                rhoedot_c=(Pedot_c/(Rgas*Te))-(Pe*Tedot_c/(Rgas*(Te**2)))
            else:
                # Edge conditions if panel is not shadowed
                Pe = P+(0.5*rho*v*v*(2-epsi)*(costheta**2))
                Pedot_Y = Pdot_Y+(0.5*rhodot_Y*v*v*(2-epsi)*(costheta**2))+(0.5*rho*v*v*(-epsidot_Y)*(costheta**2))
                Pedot_V = (rho*v*vdot_V*(2-epsi)*(costheta**2))+(0.5*rho*v*v*(-epsidot_V)*(costheta**2))
                Pedot_gam = (rho*v*v*(2-epsi)*(costheta*costhetadot_gam))
                Pedot_theta = (rho*v*v*(2-epsi)*(costheta*costhetadot_theta))
                Pedot_c = (rho*v*v*(2-epsi)*costheta*costhetadot_c)

                Te = T1*((Pe/P1)**((gamma-1)/gamma))
                Tedot_Y = (((Pe/P1)**((gamma - 1)/gamma))*T1dot_Y)+(((T1*(Pe/P1)**((gamma - 1)/gamma - 1)*(gamma - 1))/(P1*gamma))*Pedot_Y)+((-(Pe*T1*(Pe/P1)**((gamma - 1)/gamma - 1)*(gamma - 1))/(P1**2*gamma))*P1dot_Y)
                Tedot_V = (((Pe/P1)**((gamma - 1)/gamma))*T1dot_V)+(((T1*(Pe/P1)**((gamma - 1)/gamma - 1)*(gamma - 1))/(P1*gamma))*Pedot_V)+((-(Pe*T1*(Pe/P1)**((gamma - 1)/gamma - 1)*(gamma - 1))/(P1**2*gamma))*P1dot_V)
                Tedot_gam = ((T1*(Pe/P1)**((gamma - 1)/gamma - 1)*(gamma - 1))/(P1*gamma))*Pedot_gam
                Tedot_theta = ((T1*(Pe/P1)**((gamma - 1)/gamma - 1)*(gamma - 1))/(P1*gamma))*Pedot_theta
                Tedot_c = ((T1*(Pe/P1)**((gamma - 1)/gamma - 1)*(gamma - 1))/(P1*gamma))*Pedot_c

                rhoe = Pe/(Rgas*Te)
                rhoedot_Y = (Pedot_Y/(Rgas*Te))-(Pe*Tedot_Y/(Rgas*(Te**2)))
                rhoedot_V = (Pedot_V/(Rgas*Te))-(Pe*Tedot_V/(Rgas*(Te**2)))
                rhoedot_gam = (Pedot_gam/(Rgas*Te))-(Pe*Tedot_gam/(Rgas*(Te**2)))
                rhoedot_theta = (Pedot_theta/(Rgas*Te))-(Pe*Tedot_theta/(Rgas*(Te**2)))
                rhoedot_c = (Pedot_c/(Rgas*Te))-(Pe*Tedot_c/(Rgas*(Te**2)))

            # Incorporating unsteady flow using piston theory
            Vn = numpy.dot(numpy.cross(numpy.array([0,KK*w,0]),pc-CG),-n_cap) # normal velocity
            Vndot_w = numpy.dot(numpy.cross(numpy.array([0,KK*1,0]),pc-CG),-n_cap)
            Vndot_c = numpy.dot(numpy.cross(numpy.array([0,KK*w,0]),pcdot_c-CGdot_c),-n_cap) + numpy.dot(numpy.cross(numpy.array([0,KK*w,0]),pc-CG),-n_capdot_c)

            ae = sqrt(gamma*Rgas*Te) # Local speed of sound
            aedot_Y = gamma*Rgas*Tedot_Y/(2*sqrt(gamma*Rgas*Te))
            aedot_V = gamma*Rgas*Tedot_V/(2*sqrt(gamma*Rgas*Te))
            aedot_gam = gamma*Rgas*Tedot_gam/(2*sqrt(gamma*Rgas*Te))
            aedot_theta = gamma*Rgas*Tedot_theta/(2*sqrt(gamma*Rgas*Te))
            aedot_c = gamma*Rgas*Tedot_c/(2*sqrt(gamma*Rgas*Te))

            PT = Pe+(rhoe*ae*Vn) # Total local pressure
            PTdot_Y = Pedot_Y+(rhoedot_Y*ae*Vn)+(rhoe*aedot_Y*Vn)
            PTdot_V = Pedot_V+(rhoedot_V*ae*Vn)+(rhoe*aedot_V*Vn)
            PTdot_gam = Pedot_gam+(rhoedot_gam*ae*Vn)+(rhoe*aedot_gam*Vn)
            PTdot_theta = Pedot_theta+(rhoedot_theta*ae*Vn)+(rhoe*aedot_theta*Vn)
            PTdot_w = rhoe*ae*Vndot_w
            PTdot_c = Pedot_c + (rhoedot_c*ae*Vn) + (rhoe*ae*Vndot_c) + (rhoe*aedot_c*Vn)

            # Compute aerodynamic force vector acting on panel at centroid
            Force = PT*dA*n_cap
            Forcedot_Y = PTdot_Y*dA*n_cap
            Forcedot_V = PTdot_V*dA*n_cap
            Forcedot_gam = PTdot_gam*dA*n_cap
            Forcedot_theta = PTdot_theta*dA*n_cap
            Forcedot_w = PTdot_w*dA*n_cap
            Forcedot_c = (PTdot_c*dA*n_cap) + (PT*dAdot_c*n_cap) + (PT*dA*n_capdot_c)

            # Update force
            F = F+Force
            Fdot_Y = Fdot_Y+Forcedot_Y
            Fdot_V = Fdot_V+Forcedot_V
            Fdot_gam = Fdot_gam+Forcedot_gam
            Fdot_theta = Fdot_theta+Forcedot_theta
            Fdot_w = Fdot_w+Forcedot_w
            Fdot_c = Fdot_c+Forcedot_c

            # Update moment
            M = M+numpy.cross(pc-CG,Force)
            Mdot_Y = Mdot_Y+numpy.cross(pc-CG,Forcedot_Y)
            Mdot_V = Mdot_V+numpy.cross(pc-CG,Forcedot_V)
            Mdot_gam = Mdot_gam+numpy.cross(pc-CG,Forcedot_gam)
            Mdot_theta = Mdot_theta+numpy.cross(pc-CG,Forcedot_theta)
            Mdot_w = Mdot_w+numpy.cross(pc-CG,Forcedot_w)
            Mdot_c = Mdot_c+numpy.cross(pc-CG,Forcedot_c)+numpy.cross(pcdot_c-CGdot_c,Force)


    # Base of fuselage (assumed t be always shadowed)
    tau = numpy.linspace(0,1,n)

    for ii in range(n-1):
        for jj in range(n-1):
            if ii!=0:
                # Vertices of quadrilateral panels
                p1 = numpy.array([0, (a_b*b_b/sqrt((b_b*cos(thetas[jj])**2)+((a_b*sin(thetas[jj]))**2)))*tau[ii]*cos(thetas[jj]), (a_b*b_b/sqrt((b_b*cos(thetas[jj])**2)+((a_b*sin(thetas[jj]))**2)))*tau[ii]*sin(thetas[jj])])
                p2 = numpy.array([0, (a_b*b_b/sqrt((b_b*cos(thetas[jj])**2)+((a_b*sin(thetas[jj]))**2)))*tau[ii+1]*cos(thetas[jj]), (a_b*b_b/sqrt((b_b*cos(thetas[jj])**2)+((a_b*sin(thetas[jj]))**2)))*tau[ii+1]*sin(thetas[jj])])
                p3 = numpy.array([0, (a_b*b_b/sqrt((b_b*cos(thetas[jj+1])**2)+((a_b*sin(thetas[jj+1]))**2)))*tau[ii+1]*cos(thetas[jj+1]), (a_b*b_b/sqrt((b_b*cos(thetas[jj+1])**2)+((a_b*sin(thetas[jj+1]))**2)))*tau[ii+1]*sin(thetas[jj+1])])
                p4 = numpy.array([0, (a_b*b_b/sqrt((b_b*cos(thetas[jj+1])**2)+((a_b*sin(thetas[jj+1]))**2)))*tau[ii]*cos(thetas[jj+1]), (a_b*b_b/sqrt((b_b*cos(thetas[jj+1])**2)+((a_b*sin(thetas[jj+1]))**2)))*tau[ii]*sin(thetas[jj+1])])

                pc = (p1+p2+p3+p4)/4 # Centroid

                n_cap = numpy.array([1.0,0,0]) # Unit normal vector
                dA = 0.5*(numpy.linalg.norm(numpy.cross(p2-p1,p1-p3))+ \
                numpy.linalg.norm(numpy.cross(p3-p1,p1-p4))) # Area of panel
            else:
                # Vertices of triangular panels
                p1 = numpy.array([0, (a_b*b_b/sqrt((b_b*cos(thetas[jj])**2)+((a_b*sin(thetas[jj]))**2)))*tau[ii]*cos(thetas[jj]), (a_b*b_b/sqrt((b_b*cos(thetas[jj])**2)+((a_b*sin(thetas[jj]))**2)))*tau[ii]*sin(thetas[jj])])
                p2 = numpy.array([0, (a_b*b_b/sqrt((b_b*cos(thetas[jj])**2)+((a_b*sin(thetas[jj]))**2)))*tau[ii+1]*cos(thetas[jj]), (a_b*b_b/sqrt((b_b*cos(thetas[jj])**2)+((a_b*sin(thetas[jj]))**2)))*tau[ii+1]*sin(thetas[jj])])
                p3 = numpy.array([0, (a_b*b_b/sqrt((b_b*cos(thetas[jj+1])**2)+((a_b*sin(thetas[jj+1]))**2)))*tau[ii+1]*cos(thetas[jj+1]), (a_b*b_b/sqrt((b_b*cos(thetas[jj+1])**2)+((a_b*sin(thetas[jj+1]))**2)))*tau[ii+1]*sin(thetas[jj+1])])

                pc = (p1+p2+p3)/3 # Centroid

                n_cap = numpy.array([1.0,0,0]) # Unit normal vector
                dA = 0.5*numpy.linalg.norm(numpy.cross(p2-p1,p1-p3)) # Area of panel

            # Edge conditions (equal to freestream)
            Pe = P
            Pedot_Y = Pdot_Y

            Te = T
            Tedot_Y = Tdot_Y

            rhoe = Pe/(Rgas*Te)
            rhoedot_Y = (Pedot_Y/(Rgas*Te))-(Pe*Tedot_Y/(Rgas*(Te**2)))

            # Unsteady flow
            Vn = numpy.dot(numpy.cross(numpy.array([0,KK*w,0]),pc-CG),-n_cap)
            Vndot_w = numpy.dot(numpy.cross(numpy.array([0,KK*1,0]),pc-CG),-n_cap)
            Vndot_c = numpy.dot(numpy.cross(numpy.array([0,KK*w,0]),-CGdot_c),-n_cap)

            ae = sqrt(gamma*Rgas*Te)
            aedot_Y = gamma*Rgas*Tedot_Y/(2*sqrt(gamma*Rgas*Te))

            PT = Pe+(rhoe*ae*Vn)
            PTdot_Y = Pedot_Y+(rhoedot_Y*ae*Vn)+(rhoe*aedot_Y*Vn)
            PTdot_w = rhoe*ae*Vndot_w
            PTdot_c = rhoe*ae*Vndot_c

            Force = PT*dA*n_cap
            Forcedot_Y = PTdot_Y*dA*n_cap
            Forcedot_w = PTdot_w*dA*n_cap
            Forcedot_c = PTdot_c*dA*n_cap

            # Update force
            F = F+Force
            Fdot_Y = Fdot_Y+Forcedot_Y
            Fdot_w = Fdot_w+Forcedot_w
            Fdot_c = Fdot_c+Forcedot_c

            # Update moment
            M = M+numpy.cross(pc-CG,Force)
            Mdot_Y = Mdot_Y+numpy.cross(pc-CG,Forcedot_Y)
            Mdot_w = Mdot_w+numpy.cross(pc-CG,Forcedot_w)
            Mdot_c = Mdot_c+numpy.cross(pc-CG,Forcedot_c)+numpy.cross(-CGdot_c,Force)

    # Control surface
    yDim = 0.3 # Length (determines span of elevons)
    aDim = 0.2 # Semi-major axis of elliptical airfoil (determines chord length)
    bDim = 0.2/4 # Semi-minor axis of elliptical airfoil

    thetas = numpy.linspace(0,2*pi,n)
    ys = numpy.linspace(a_b, a_b+yDim, n)

    for idx1 in range(n-1):
        for idx2 in range(n-1):
            # Vertices of panels
            R1 = numpy.array([aDim*cos(thetas[idx1]), ys[idx2], bDim*sin(thetas[idx1])])
            R2 = numpy.array([aDim*cos(thetas[idx1+1]), ys[idx2], bDim*sin(thetas[idx1+1])])
            R3 = numpy.array([aDim*cos(thetas[idx1+1]), ys[idx2+1], bDim*sin(thetas[idx1+1])])
            R4 = numpy.array([aDim*cos(thetas[idx1]), ys[idx2+1], bDim*sin(thetas[idx1])])

            RcBody = (R1+R2+R3+R4)/4
            dA = 0.5*(numpy.linalg.norm(numpy.cross(R2-R1,R1-R3))+numpy.linalg.norm(numpy.cross(R3-R1,R1-R4)))

            pc = (RcBody[0]*numpy.array([cos(delta), 0, -sin(delta)])) + (RcBody[1]*numpy.array([0, 1.0, 0])) \
            + (RcBody[2]*numpy.array([sin(delta), 0, cos(delta)])) # Centroid

            pcdot_w = (RcBody[0]*numpy.array([-sin(delta),0,-cos(delta)])*deldot_w) + \
            (RcBody[2]*numpy.array([cos(delta),0,-sin(delta)])*deldot_w)

            nBar = numpy.cross(R2-R1, R3-R1)
            nCapBody = nBar/numpy.linalg.norm(nBar)
            n_cap = (nCapBody[0]*numpy.array([cos(delta),0,-sin(delta)])) + (nCapBody[1]*numpy.array([0,1.0,0])) \
            + (nCapBody[2]*numpy.array([sin(delta),0,cos(delta)])) # Unit normal

            n_capdot_w = (nCapBody[0]*numpy.array([-sin(delta),0,-cos(delta)])*deldot_w) \
            + (nCapBody[2]*numpy.array([cos(delta),0,-sin(delta)])*deldot_w)

            # Angle beween inward-pointing unit normal and relative wind
            costheta = numpy.dot(W,n_cap)
            costhetadot_theta = numpy.dot(Wdot_theta,n_cap)
            costhetadot_gam = numpy.dot(Wdot_gam,n_cap)
            costhetadot_w = numpy.dot(W,n_capdot_w)

            if costheta<0:
                # If shadowed
                Pe = P
                Pedot_Y = Pdot_Y
                Pedot_V = 0
                Pedot_gam = 0
                Pedot_theta = 0
                Pedot_w=0

                Te = T
                Tedot_Y = Tdot_Y
                Tedot_V = 0
                Tedot_gam = 0
                Tedot_theta = 0
                Tedot_w = 0

                rhoe = Pe/(Rgas*Te)
                rhoedot_Y = (Pedot_Y/(Rgas*Te))-(Pe*Tedot_Y/(Rgas*(Te**2)))
                rhoedot_V = (Pedot_V/(Rgas*Te))-(Pe*Tedot_V/(Rgas*(Te**2)))
                rhoedot_gam = (Pedot_gam/(Rgas*Te))-(Pe*Tedot_gam/(Rgas*(Te**2)))
                rhoedot_theta = (Pedot_theta/(Rgas*Te))-(Pe*Tedot_theta/(Rgas*(Te**2)))
                rhoedot_w = (Pedot_w/(Rgas*Te))-(Pe*Tedot_w/(Rgas*(Te**2)))
            else:
                # If not shadowed
                Pe = P+(0.5*rho*v*v*(2-epsi)*(costheta**2))
                Pedot_Y = Pdot_Y+(0.5*rhodot_Y*v*v*(2-epsi)*(costheta**2))+(0.5*rho*v*v*(-epsidot_Y)*(costheta**2))
                Pedot_V = (rho*v*vdot_V*(2-epsi)*(costheta**2))+(0.5*rho*v*v*(-epsidot_V)*(costheta**2))
                Pedot_gam = (rho*v*v*(2-epsi)*(costheta*costhetadot_gam))
                Pedot_theta = (rho*v*v*(2-epsi)*(costheta*costhetadot_theta))
                Pedot_w = (0.5*rho*v*v*(2-epsi)*(2*costheta*costhetadot_w))

                Te = T1*((Pe/P1)**((gamma-1)/gamma))
                Tedot_Y = (((Pe/P1)**((gamma - 1)/gamma))*T1dot_Y)+(((T1*(Pe/P1)**((gamma - 1)/gamma - 1)*(gamma - 1))/(P1*gamma))*Pedot_Y)+((-(Pe*T1*(Pe/P1)**((gamma - 1)/gamma - 1)*(gamma - 1))/(P1**2*gamma))*P1dot_Y)
                Tedot_V = (((Pe/P1)**((gamma - 1)/gamma))*T1dot_V)+(((T1*(Pe/P1)**((gamma - 1)/gamma - 1)*(gamma - 1))/(P1*gamma))*Pedot_V)+((-(Pe*T1*(Pe/P1)**((gamma - 1)/gamma - 1)*(gamma - 1))/(P1**2*gamma))*P1dot_V)
                Tedot_gam = ((T1*(Pe/P1)**((gamma - 1)/gamma - 1)*(gamma - 1))/(P1*gamma))*Pedot_gam
                Tedot_theta = ((T1*(Pe/P1)**((gamma - 1)/gamma - 1)*(gamma - 1))/(P1*gamma))*Pedot_theta
                Tedot_w = ((T1*(Pe/P1)**((gamma - 1)/gamma - 1)*(gamma - 1))/(P1*gamma))*Pedot_w

                rhoe = Pe/(Rgas*Te)
                rhoedot_Y = (Pedot_Y/(Rgas*Te))-(Pe*Tedot_Y/(Rgas*(Te**2)))
                rhoedot_V = (Pedot_V/(Rgas*Te))-(Pe*Tedot_V/(Rgas*(Te**2)))
                rhoedot_gam = (Pedot_gam/(Rgas*Te))-(Pe*Tedot_gam/(Rgas*(Te**2)))
                rhoedot_theta = (Pedot_theta/(Rgas*Te))-(Pe*Tedot_theta/(Rgas*(Te**2)))
                rhoedot_w = (Pedot_w/(Rgas*Te))-(Pe*Tedot_w/(Rgas*(Te**2)))

            # Accounting for unsteady flow using piston theory
            Vn = numpy.dot(numpy.cross(numpy.array([0,KK*w,0]),pc-CG),-n_cap)
            Vndot_w = numpy.dot(numpy.cross(numpy.array([0,KK*1,0]),pc-CG),-n_cap) + numpy.dot(numpy.cross(numpy.array([0,KK*w,0]),pc-CG),-n_capdot_w) + numpy.dot(numpy.cross(numpy.array([0,KK*w,0]),pcdot_w),-n_cap)
            Vndot_c = numpy.dot(numpy.cross(numpy.array([0,KK*w,0]),-CGdot_c),-n_cap)

            ae = sqrt(gamma*Rgas*Te)
            aedot_Y = gamma*Rgas*Tedot_Y/(2*sqrt(gamma*Rgas*Te))
            aedot_V = gamma*Rgas*Tedot_V/(2*sqrt(gamma*Rgas*Te))
            aedot_gam = gamma*Rgas*Tedot_gam/(2*sqrt(gamma*Rgas*Te))
            aedot_theta = gamma*Rgas*Tedot_theta/(2*sqrt(gamma*Rgas*Te))
            aedot_w = gamma*Rgas*Tedot_w/(2*sqrt(gamma*Rgas*Te))

            PT = Pe+(rhoe*ae*Vn)
            PTdot_Y = Pedot_Y+(rhoedot_Y*ae*Vn)+(rhoe*aedot_Y*Vn)
            PTdot_V = Pedot_V+(rhoedot_V*ae*Vn)+(rhoe*aedot_V*Vn)
            PTdot_gam = Pedot_gam+(rhoedot_gam*ae*Vn)+(rhoe*aedot_gam*Vn)
            PTdot_theta = Pedot_theta+(rhoedot_theta*ae*Vn)+(rhoe*aedot_theta*Vn)
            PTdot_w = (rhoe*ae*Vndot_w) + Pedot_w + (rhoedot_w*ae*Vn)+(rhoe*aedot_w*Vn)
            PTdot_c = rhoe*ae*Vndot_c

            Force = PT*dA*n_cap
            Forcedot_Y = PTdot_Y*dA*n_cap
            Forcedot_V = PTdot_V*dA*n_cap
            Forcedot_gam = PTdot_gam*dA*n_cap
            Forcedot_theta = PTdot_theta*dA*n_cap
            Forcedot_w = (PTdot_w*dA*n_cap) + (PT*dA*n_capdot_w)
            Forcedot_c = PTdot_c*dA*n_cap

            # Update force
            F = F + Force
            Fdot_Y = Fdot_Y + Forcedot_Y
            Fdot_V = Fdot_V + Forcedot_V
            Fdot_gam = Fdot_gam + Forcedot_gam
            Fdot_theta = Fdot_theta + Forcedot_theta
            Fdot_w = Fdot_w + Forcedot_w
            Fdot_c = Fdot_c + Forcedot_c

            # Update moment
            M = M + numpy.cross(pc-CG,Force)
            Mdot_Y = Mdot_Y + numpy.cross(pc-CG,Forcedot_Y)
            Mdot_V = Mdot_V + numpy.cross(pc-CG,Forcedot_V)
            Mdot_gam = Mdot_gam + numpy.cross(pc-CG,Forcedot_gam)
            Mdot_theta = Mdot_theta + numpy.cross(pc-CG,Forcedot_theta)
            Mdot_w = Mdot_w + numpy.cross(pc-CG,Forcedot_w) + numpy.cross(pcdot_w, Force)
            Mdot_c = Mdot_c + numpy.cross(pc-CG,Forcedot_c) + numpy.cross(-CGdot_c,Force)

    # Extract desired force and moment components and their derivatives.
    # Can be updated to include all components if using 6-DOF dynamics.
    fx = 2*F[0]
    fz = 2*F[2]
    My = 2*M[1]

    fxdot_Y = 2*Fdot_Y[0]
    fzdot_Y = 2*Fdot_Y[2]
    Mydot_Y = 2*Mdot_Y[1]

    fxdot_V = 2*Fdot_V[0]
    fzdot_V = 2*Fdot_V[2]
    Mydot_V = 2*Mdot_V[1]

    fxdot_gam = 2*Fdot_gam[0]
    fzdot_gam = 2*Fdot_gam[2]
    Mydot_gam = 2*Mdot_gam[1]

    fxdot_theta = 2*Fdot_theta[0]
    fzdot_theta = 2*Fdot_theta[2]
    Mydot_theta = 2*Mdot_theta[1]

    fxdot_w = 2*Fdot_w[0]
    fzdot_w = 2*Fdot_w[2]
    Mydot_w = 2*Mdot_w[1]

    fxdot_c = 2*Fdot_c[0]
    fzdot_c = 2*Fdot_c[2]
    Mydot_c = 2*Mdot_c[1]

    # Construct the return vector
    aerocoeff = numpy.array([fx,fz,My,fxdot_Y,fzdot_Y,Mydot_Y,fxdot_V,fzdot_V,Mydot_V,fxdot_gam,fzdot_gam,Mydot_gam,fxdot_theta,fzdot_theta,Mydot_theta,fxdot_w,fzdot_w,Mydot_w,fxdot_c,fzdot_c,Mydot_c])
    return aerocoeff

def aero_autodiff_u(theta,gam,V,Y,w,c,u_sat):
    '''
    This function calculates the derivatives of the desired aerodynamic force
    components with respect to the fictitious control, using automatic
    differentiation. These values are returned in a numpy array.

    The vehicle of interest possesses an elliptic paraboloid geometry.

    The aerodynamic forces and moments are modeled using panel methods and
    Modified Newtonian flow theory. Unsteady flow is modeled using piston theory.

    Inputs:
    theta = downrange angle (longitude), rad
    gam = flight-path-angle, rad
    V = scaled velocity, nd
    Y = scaled altitude, nd
    w = pitch rate, rad/s
    c = length of fuselage, m
    u_sat = fictitious control used for bounding original control,
            elevon deflection command (del_c)
    '''

    # Angle-of-attack
    alpha = theta-gam

    # Dimensionalize velocity
    v0 = 3000.0 # Velocity scale factor
    y0 = 40000.0 # Altitude scale factor
    v = V*v0

    # Compute elevon deflection angle from fictitious control
    ub_del_c = 90*pi/180
    lb_del_c = -ub_del_c
    slp_u = 0.3

    del_c = (ub_del_c-((ub_del_c-lb_del_c)/(1+exp(slp_u*u_sat))))
    del_c_dot_u_sat = (-(slp_u*exp(slp_u*u_sat)*(lb_del_c - ub_del_c))/(exp(slp_u*u_sat) + 1)**2)

    k=1
    KK=1

    delta = del_c+(k*w)
    deldot_u_sat = del_c_dot_u_sat

    # Dimensions of fuselage base
    a_b = 0.5 # Semi-major axis
    b_b = 0.5 # Semi-minor axis

    # Fuselage
    rBaseFuselage = 0.5
    densityFuselage = 800

    # Warhead
    densityWarhead = 17000
    rWarhead = 0.2
    warheadPosition = 2.5

    # Subsystem
    lengthSubsystem = 2.2
    rBaseSubsystem = 0.4
    densitySubsystem = 100

    # Mass
    massFuselage = (densityFuselage*pi*(rBaseFuselage**2)*c/2) - \
    (densityFuselage*(4/3)*pi*(rWarhead**3)) - (densityFuselage*pi*(rBaseSubsystem**2)*lengthSubsystem/2)

    massWarhead = densityWarhead*(4/3)*pi*(rWarhead**3)
    massSubsystem = (densitySubsystem*pi*(rBaseSubsystem**2)*lengthSubsystem/2)
    mass = massFuselage + massWarhead + massSubsystem

    # Center of mass
    cgFuselage = (((densityFuselage*pi*(rBaseFuselage**2)*c/2)*c/3) \
    - (densityFuselage*(4/3)*pi*(rWarhead**3)*warheadPosition) \
    - ((densityFuselage*pi*(rBaseSubsystem**2)*lengthSubsystem/2)*lengthSubsystem/3))/massFuselage

    cgWarhead = warheadPosition
    cgSubsystem = lengthSubsystem/3
    cg = ((massFuselage*cgFuselage) + (massWarhead*cgWarhead) + (massSubsystem*cgSubsystem))/mass

    CG = numpy.array([cg,0,0])

    # Freestream conditions

    Rgas = 287
    gamma = 1.4
    rho0 = 1.225
    HScale = 8500

    W = numpy.array([-cos(alpha),0,-sin(alpha)])

    rho = rho0*exp(-(Y*y0/HScale))
    P = 101325*exp(-(Y*y0/HScale))
    T = P/(rho*Rgas)

    Mach = v/sqrt(gamma*Rgas*T)
    epsi=(((gamma-1)*(Mach**2))+2)/((gamma+1)*(Mach**2))

    # Downstream conditions behind normal portion of shock
    P1 = P*(((2*gamma*(Mach**2))-(gamma-1))/(gamma+1))
    rho1 = rho/epsi
    T1 = P1/(rho1*Rgas)

    # Initialize force and moment vectors
    F = numpy.array([0,0,0])
    Fdot_u_sat = numpy.array([0,0,0])

    M = numpy.array([0,0,0])
    Mdot_u_sat = numpy.array([0,0,0])

    # Control surface
    yDim = 0.3
    aDim = 0.2
    bDim = 0.2/4

    n = 11
    thetas = numpy.linspace(0,2*pi,n)
    ys = numpy.linspace(a_b, a_b+yDim, n)

    for idx1 in range(n-1):
        for idx2 in range(n-1):
            # Calculate vertices of panels
            R1 = numpy.array([aDim*cos(thetas[idx1]),ys[idx2],bDim*sin(thetas[idx1])])
            R2 = numpy.array([aDim*cos(thetas[idx1+1]),ys[idx2],bDim*sin(thetas[idx1+1])])
            R3 = numpy.array([aDim*cos(thetas[idx1+1]),ys[idx2+1],bDim*sin(thetas[idx1+1])])
            R4 = numpy.array([aDim*cos(thetas[idx1]),ys[idx2+1],bDim*sin(thetas[idx1])])

            RcBody = (R1+R2+R3+R4)/4 # Centroid of panel
            dA = 0.5*(numpy.linalg.norm(numpy.cross(R2-R1,R1-R3))+numpy.linalg.norm(numpy.cross(R3-R1,R1-R4))) # Area of panel

            # Centroid of panel in body frame
            pc = (RcBody[0]*numpy.array([cos(delta),0,-sin(delta)])) + (RcBody[1]*numpy.array([0,1.0,0])) \
            + (RcBody[2]*numpy.array([sin(delta),0,cos(delta)]))

            pcdot_u_sat = (RcBody[0]*numpy.array([-sin(delta),0,-cos(delta)])*deldot_u_sat) + \
            (RcBody[2]*numpy.array([cos(delta),0,-sin(delta)])*deldot_u_sat)

            # Normal vector computation
            nBar = numpy.cross(R2-R1, R3-R1)
            nCapBody = nBar/numpy.linalg.norm(nBar)

            n_cap = (nCapBody[0]*numpy.array([cos(delta),0,-sin(delta)])) + (nCapBody[1]*numpy.array([0,1.0,0])) \
            + (nCapBody[2]*numpy.array([sin(delta),0,cos(delta)]))

            n_capdot_u_sat = (nCapBody[0]*numpy.array([-sin(delta),0,-cos(delta)])*deldot_u_sat) \
            + (nCapBody[2]*numpy.array([cos(delta),0,-sin(delta)])*deldot_u_sat)

            # Angle between inward-pointing unit normal and relative wind
            costheta = numpy.dot(W,n_cap)
            costhetadot_u_sat = numpy.dot(W,n_capdot_u_sat)

            if costheta<0:
                # Edge conditions if shadowed
                Pe = P
                Pedot_u_sat = 0

                Te = T
                Tedot_u_sat = 0

                rhoe = Pe/(Rgas*Te)
                rhoedot_u_sat = (Pedot_u_sat/(Rgas*Te))-(Pe*Tedot_u_sat/(Rgas*(Te**2)))
            else:
                # Edge conditions if not shadowed
                Pe = P+(0.5*rho*v*v*(2-epsi)*(costheta**2))
                Pedot_u_sat = (0.5*rho*v*v*(2-epsi)*(2*costheta*costhetadot_u_sat))

                Te = T1*((Pe/P1)**((gamma-1)/gamma))
                Tedot_u_sat = ((T1*(Pe/P1)**((gamma - 1)/gamma - 1)*(gamma - 1))/(P1*gamma))*Pedot_u_sat

                rhoe = Pe/(Rgas*Te)
                rhoedot_u_sat = (Pedot_u_sat/(Rgas*Te))-(Pe*Tedot_u_sat/(Rgas*(Te**2)))

            # Accounting for unsteady flow using piston theory
            Vn = numpy.dot(numpy.cross(numpy.array([0,KK*w,0]),pc-CG),-n_cap) # Normal velocity
            Vndot_u_sat = numpy.dot(numpy.cross(numpy.array([0,KK*w,0]),pc-CG),-n_capdot_u_sat) + numpy.dot(numpy.cross(numpy.array([0,KK*w,0]),pcdot_u_sat),-n_cap)

            ae = sqrt(gamma*Rgas*Te) # Local speed of sound
            aedot_u_sat = gamma*Rgas*Tedot_u_sat/(2*sqrt(gamma*Rgas*Te))

            PT = Pe+(rhoe*ae*Vn) # Total pressure
            PTdot_u_sat = (rhoe*ae*Vndot_u_sat) + Pedot_u_sat + (rhoedot_u_sat*ae*Vn)+(rhoe*aedot_u_sat*Vn)

            Force = PT*dA*n_cap
            Forcedot_u_sat = (PTdot_u_sat*dA*n_cap) + (PT*dA*n_capdot_u_sat)

            # Update force and moment
            F = F + Force;
            Fdot_u_sat = Fdot_u_sat + Forcedot_u_sat

            M = M + numpy.cross(pc-CG,Force)
            Mdot_u_sat = Mdot_u_sat + numpy.cross(pc-CG,Forcedot_u_sat) + numpy.cross(pcdot_u_sat, Force)

    # Extract desired quantities
    fxdot_u_sat = 2*Fdot_u_sat[0]
    fzdot_u_sat = 2*Fdot_u_sat[2]
    Mydot_u_sat = 2*Mdot_u_sat[1]

    # COnstruct the return vector
    aerocoeff = numpy.array([fxdot_u_sat,fzdot_u_sat,Mydot_u_sat])
    return aerocoeff

states = numpy.array([0.5,0,0.8,0,0,0,4,0.1,0.1,0.1,0.1,0.1,0.1])
params = numpy.array([100,0.1])
tau = 0.7
hCont = 0
sCont = 0
flag = 1
epsi = 0.1

#XDots = odes(tau,states,params,hCont, sCont, epsi, flag)
#print(XDots)

#print(fsolve(controlFun, 0.1, args=(states, tau*100, epsi), fprime=None, full_output=0, \
#col_deriv=0, xtol=1e-3, maxfev=0, band=None, epsfcn=None, factor=100, diag=None))

#print(controlFun(0.1,states, tau*100, epsi))
#print(aero_autodiff_u(0,0,0.5,0.8,0,4,0.1))

#bc = bcs(states,states,params,hCont, sCont, epsi, flag)
#print(bc)
