from math import *
import numpy
from scipy.optimize import bisect

def odes(tau, states, params, aux):

    v = states[0]
    gam = states[1]
    h = states[2]
    s = states[3]
    c = states[4]

    lamv = states[5]
    lamgam = states[6]
    lamh = states[7]
    lams = states[8]
    lamc = states[9]

    tf = states[10]

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

    # Subsystem
    lengthSubsystem = 2.2
    rBaseSubsystem = 0.4
    densitySubsystem = 100.0

    alpha = bisect(controlFun, -(pi/3), pi/3, args=(states, t), xtol=1e-3, rtol=1e-3, maxiter=10000, full_output=False, disp=False)

    aeroCoeff = aero_autodiff(v,h,c,alpha)

    fx = aeroCoeff[0]
    fz = aeroCoeff[1]

    fxdot_h = aeroCoeff[2]
    fzdot_h = aeroCoeff[3]

    fxdot_v = aeroCoeff[4]
    fzdot_v = aeroCoeff[5]

    fxdot_c = aeroCoeff[6]
    fzdot_c = aeroCoeff[7]

    vdot = ((cos(alpha)*fx + sin(alpha)*fz)/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2) - (mu*sin(gam))/(R + h*h0)**2)/v0;
    gamdot = (v*v0*cos(gam))/(R + h*h0) - (cos(alpha)*fz - sin(alpha)*fx)/(v*v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) - (mu*cos(gam))/(v*v0*(R + h*h0)**2);
    hdot = (v*v0*sin(gam))/h0;
    sdot = (v*v0*cos(gam))/(R + h*h0);
    cdot = -(2**(1/2)*k*v**3*v0**3*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2))/(QStar*densityFuselage);

    lamvdot = (3*2**(1/2)*k*lamc*v**2*v0**3*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2))/(QStar*densityFuselage) - (lams*v0*cos(gam))/(R + h*h0) - (lamv*(cos(alpha)*fxdot_v + sin(alpha)*fzdot_v))/(v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) - (lamh*v0*sin(gam))/h0 - lamgam*((v0*cos(gam))/(R + h*h0) - (cos(alpha)*fzdot_v - sin(alpha)*fxdot_v)/(v*v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) + (cos(alpha)*fz - sin(alpha)*fx)/(v**2*v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) + (mu*cos(gam))/(v**2*v0*(R + h*h0)**2))

    lamgamdot = lamgam*((v*v0*sin(gam))/(R + h*h0) - (mu*sin(gam))/(v*v0*(R + h*h0)**2)) + (lamv*mu*cos(gam))/(v0*(R + h*h0)**2) - (lamh*v*v0*cos(gam))/h0 + (lams*v*v0*sin(gam))/(R + h*h0)

    lamhdot = (h0*lams*v*v0*cos(gam))/(R + h*h0)**2 - (lamv*((sin(alpha)*fzdot_h + cos(alpha)*fxdot_h)/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2) + (2*h0*mu*sin(gam))/(R + h*h0)**3))/v0 - lamgam*((sin(alpha)*fxdot_h - cos(alpha)*fzdot_h)/(v*v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) - (h0*v*v0*cos(gam))/(R + h*h0)**2 + (2*h0*mu*cos(gam))/(v*v0*(R + h*h0)**3)) - (2**(1/2)*c*h0*k*lamc*rho0*v**3*v0**3*exp(-(h*h0)/HScale))/(2*HScale*QStar*densityFuselage*rBaseFuselage**2*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2))

    lamsdot = 0

    lamcdot = (2**(1/2)*k*lamc*rho0*v**3*v0**3*exp(-(h*h0)/HScale))/(2*QStar*densityFuselage*rBaseFuselage**2*((c*rho0*exp(-(h*h0)/HScale))/rBaseFuselage**2)**(1/2)) - (lamv*((sin(alpha)*fzdot_c + cos(alpha)*fxdot_c)/((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2) - (densityFuselage*rBaseFuselage**2*pi*(cos(alpha)*fx + sin(alpha)*fz))/(2*((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2)))/v0 - lamgam*((sin(alpha)*fxdot_c - cos(alpha)*fzdot_c)/(v*v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) + (densityFuselage*rBaseFuselage**2*pi*(cos(alpha)*fz - sin(alpha)*fx))/(2*v*v0*((4*densityWarhead*pi*rWarhead**3)/3 - (4*densityFuselage*pi*rWarhead**3)/3 + (c*densityFuselage*pi*rBaseFuselage**2)/2 - (densityFuselage*lengthSubsystem*pi*rBaseSubsystem**2)/2 + (densitySubsystem*lengthSubsystem*pi*rBaseSubsystem**2)/2)**2))

    tfdot = 0
    XDot = tf*numpy.array([vdot,gamdot,hdot,sdot,cdot,lamvdot,lamgamdot,lamhdot,lamsdot,lamcdot,tfdot])
    return XDot

def controlFun(alpha, states, t):
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

    # Subsystem
    lengthSubsystem = 2.2
    rBaseSubsystem = 0.4
    densitySubsystem = 100.0

    v = states[0]
    gam = states[1]
    h = states[2]
    s = states[3]
    c = states[4]

    lamv = states[5]
    lamgam = states[6]
    lamh = states[7]
    lams = states[8]
    lamc = states[9]

    aeroCoeff = aero_autodiff_u(v,h,c,alpha)
    fx = aeroCoeff[0]
    fz = aeroCoeff[1]
    fxdot_alpha = aeroCoeff[2]
    fzdot_alpha = aeroCoeff[3]

    dHdu = (lamv*(sin(alpha)*fzdot_alpha + cos(alpha)*fz - sin(alpha)*fx + cos(alpha)*fxdot_alpha))/(v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2)) + (lamgam*(sin(alpha)*fxdot_alpha + cos(alpha)*fx + sin(alpha)*fz - cos(alpha)*fzdot_alpha))/(v*v0*((4*densityWarhead*rWarhead**3*pi)/3 - (4*densityFuselage*rWarhead**3*pi)/3 + (c*densityFuselage*rBaseFuselage**2*pi)/2 - (densityFuselage*lengthSubsystem*rBaseSubsystem**2*pi)/2 + (densitySubsystem*lengthSubsystem*rBaseSubsystem**2*pi)/2));

    return dHdu

def bcs(Ya, Yb, params, aux):
    hCont = aux[0]
    sCont = aux[1]
    flag = aux[2]

    v0 = Ya[0]
    gam0 = Ya[1]
    h0 = Ya[2]
    s0 = Ya[3]
    c0 = Ya[4]
    lamv0 = Ya[5]
    lamgam0 = Ya[6]
    lamh0 = Ya[7]
    lams0 = Ya[8]
    lamc0 = Ya[9]
    lamtf0 = Ya[10]

    v = Yb[0]
    gam = Yb[1]
    h = Yb[2]
    s = Yb[3]
    c = Yb[4]
    lamv = Yb[5]
    lamgam = Yb[6]
    lamh = Yb[7]
    lams = Yb[8]
    lamc = Yb[9]
    tf = Yb[10]

    Xdot = odes(1, Yb, params, aux)

    weight = 0

    Ham = (weight*1) + (numpy.dot(numpy.array([lamv,lamgam,lamh,lams,lamc]),Xdot[:5])/tf)

    if flag == 1:
        zeroVec = numpy.array([v0 - 1,\
        gam0 - 0,\
        h0 - 1,\
        s0 - 0,\
        c0 - 4,\
        lamv + v,\
        lamgam - 0,\
        lamh - 0,\
        s - sCont,\
        lamc - 0,\
        Ham - 0])
    else:
        zeroVec = numpy.array([v0 - 1,\
        gam0 - 0,\
        h0 - 1,\
        s0 - 0,\
        c0 - 4,\
        lamv + v,\
        lamgam - 0,\
        h - hCont,\
        s - sCont,\
        lamc - 0,\
        Ham - 0])

    return zeroVec

def aero_autodiff(V,Y,c,alpha):
    a_b=0.5
    b_b=0.5

    n=11

    v0 = 3000.0
    y0 = 40000.0

    v=V*v0
    vdot_V=v0

    Rgas=287.0
    rho0 = 1.225
    HScale = 8500.0
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

    W=numpy.array([-cos(alpha),0,-sin(alpha)])

    F=numpy.array([0,0,0])
    Fdot_V=numpy.array([0,0,0])
    Fdot_Y=numpy.array([0,0,0])
    Fdot_c=numpy.array([0,0,0])

    x=numpy.linspace(0,c,n)
    xdot_c = numpy.linspace(0,1,n)
    thetas=numpy.linspace(pi/2,3*pi/2,n)

    for ii in range(n-1):
        for jj in range(n-1):
            if ii!=n-2:
                X=x[ii]
                Xdot_c = xdot_c[ii]

                Xp=x[ii+1]
                Xpdot_c = xdot_c[ii+1]

                thetass=thetas[jj]
                thetap=thetas[jj+1]

                p1=numpy.array([X,a_b*sqrt(1-(X/c))*cos(thetass),b_b*sqrt(1-(X/c))*sin(thetass)])
                p1dot_c = numpy.array([Xdot_c,a_b*(-(Xdot_c/c)+(X/(c**2)))*cos(thetass)/(2*sqrt(1-(X/c))),b_b*(-(Xdot_c/c)+(X/(c**2)))*sin(thetass)/(2*sqrt(1-(X/c)))])

                p2=numpy.array([Xp,a_b*sqrt(1-(Xp/c))*cos(thetass),b_b*sqrt(1-(Xp/c))*sin(thetass)])
                p2dot_c = numpy.array([Xpdot_c,a_b*(-(Xpdot_c/c)+(Xp/(c**2)))*cos(thetass)/(2*sqrt(1-(Xp/c))),b_b*(-(Xpdot_c/c)+(Xp/(c**2)))*sin(thetass)/(2*sqrt(1-(Xp/c)))])

                p3=numpy.array([Xp,a_b*sqrt(1-(Xp/c))*cos(thetap),b_b*sqrt(1-(Xp/c))*sin(thetap)])
                p3dot_c = numpy.array([Xpdot_c,a_b*(-(Xpdot_c/c)+(Xp/(c**2)))*cos(thetap)/(2*sqrt(1-(Xp/c))),b_b*(-(Xpdot_c/c)+(Xp/(c**2)))*sin(thetap)/(2*sqrt(1-(Xp/c)))])

                p4=numpy.array([X,a_b*sqrt(1-(X/c))*cos(thetap),b_b*sqrt(1-(X/c))*sin(thetap)])
                p4dot_c = numpy.array([Xdot_c,a_b*(-(Xdot_c/c)+(X/(c**2)))*cos(thetap)/(2*sqrt(1-(X/c))),b_b*(-(Xdot_c/c)+(X/(c**2)))*sin(thetap)/(2*sqrt(1-(X/c)))])

                N=numpy.cross(p2-p1,p3-p1)
                Ndot_c = numpy.cross(p2dot_c-p1dot_c,p3-p1) + numpy.cross(p2-p1,p3dot_c-p1dot_c)

                normN = numpy.linalg.norm(N)
                normNdot_c = ((N[0]*Ndot_c[0]) + (N[1]*Ndot_c[1]) + (N[2]*Ndot_c[2]))/normN

                n_cap=N/numpy.linalg.norm(N)
                n_capdot_c = (Ndot_c/normN) - (N*normNdot_c/(normN**2))

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
                X=x[ii]
                Xdot_c = xdot_c[ii]

                Xp=x[ii+1]
                Xpdot_c = xdot_c[ii+1]

                thetass=thetas[jj]
                thetap=thetas[jj+1]

                p1=numpy.array([X,a_b*sqrt(1-(X/c))*cos(thetass),b_b*sqrt(1-(X/c))*sin(thetass)])
                p1dot_c = numpy.array([Xdot_c,a_b*(-(Xdot_c/c)+(X/(c**2)))*cos(thetass)/(2*sqrt(1-(X/c))),b_b*(-(Xdot_c/c)+(X/(c**2)))*sin(thetass)/(2*sqrt(1-(X/c)))])

                p2=numpy.array([Xp,a_b*sqrt(1-(Xp/c))*cos(thetass),b_b*sqrt(1-(Xp/c))*sin(thetass)])

                if ii == n-2:
                    p2dot_c = numpy.array([Xpdot_c,0,0])
                else:
                    p2dot_c = numpy.array([Xpdot_c,a_b*(-(Xpdot_c/c)+(Xp/(c**2)))*cos(thetass)/(2*sqrt(1-(Xp/c))),b_b*(-(Xpdot_c/c)+(Xp/(c**2)))*sin(thetass)/(2*sqrt(1-(Xp/c)))])

                p3=numpy.array([X,a_b*sqrt(1-(X/c))*cos(thetap),b_b*sqrt(1-(X/c))*sin(thetap)])
                p3dot_c = numpy.array([Xdot_c,a_b*(-(Xdot_c/c)+(X/(c**2)))*cos(thetap)/(2*sqrt(1-(X/c))),b_b*(-(Xdot_c/c)+(X/(c**2)))*sin(thetap)/(2*sqrt(1-(X/c)))])

                N=numpy.cross(p2-p1,p3-p1)
                Ndot_c = numpy.cross(p2dot_c-p1dot_c,p3-p1) + numpy.cross(p2-p1,p3dot_c-p1dot_c)

                normN = numpy.linalg.norm(N)
                normNdot_c = ((N[0]*Ndot_c[0]) + (N[1]*Ndot_c[1]) + (N[2]*Ndot_c[2]))/normN

                n_cap=N/normN
                n_capdot_c = (Ndot_c/normN) - (N*normNdot_c/(normN**2))

                cross1 = numpy.cross(p2-p1,p1-p3)
                cross1dot_c = numpy.cross(p2dot_c-p1dot_c,p1-p3) + numpy.cross(p2-p1,p1dot_c-p3dot_c)

                normcross1 = numpy.linalg.norm(cross1)
                normcross1dot_c = ((cross1[0]*cross1dot_c[0]) + (cross1[1]*cross1dot_c[1]) + (cross1[2]*cross1dot_c[2]))/normcross1

                dA=0.5*normcross1
                dAdot_c = 0.5*normcross1dot_c


            costheta=numpy.dot(W,n_cap)
            costhetadot_c = numpy.dot(W,n_capdot_c)

            if costheta<0:
                Pe=P
                Pedot_Y=Pdot_Y
                Pedot_V=0
                Pedot_c = 0
            else:
                Pe=P+(0.5*rho*v*v*(2-epsi)*(costheta**2))
                Pedot_Y=Pdot_Y+(0.5*rhodot_Y*v*v*(2-epsi)*(costheta**2))+(0.5*rho*v*v*(-epsidot_Y)*(costheta**2))
                Pedot_V=(rho*v*vdot_V*(2-epsi)*(costheta**2))+(0.5*rho*v*v*(-epsidot_V)*(costheta**2))
                Pedot_c=(rho*v*v*(2-epsi)*costheta*costhetadot_c)

            PT=Pe
            PTdot_Y=Pedot_Y
            PTdot_V=Pedot_V
            PTdot_c=Pedot_c

            Force=PT*dA*n_cap
            Forcedot_Y=PTdot_Y*dA*n_cap
            Forcedot_V=PTdot_V*dA*n_cap
            Forcedot_c=(PTdot_c*dA*n_cap) + (PT*dAdot_c*n_cap) + (PT*dA*n_capdot_c)

            F=F+Force
            Fdot_Y=Fdot_Y+Forcedot_Y
            Fdot_V=Fdot_V+Forcedot_V
            Fdot_c=Fdot_c+Forcedot_c

    #Base
    n=11
    tau=numpy.linspace(0,1,n)

    for ii in range(n-1):
        for jj in range(n-1):
            if ii!=0:
                p1=numpy.array([0,(a_b*b_b/sqrt((b_b*cos(thetas[jj])**2)+((a_b*sin(thetas[jj]))**2)))*tau[ii]*cos(thetas[jj]),(a_b*b_b/sqrt((b_b*cos(thetas[jj])**2)+((a_b*sin(thetas[jj]))**2)))*tau[ii]*sin(thetas[jj])])
                p2=numpy.array([0,(a_b*b_b/sqrt((b_b*cos(thetas[jj])**2)+((a_b*sin(thetas[jj]))**2)))*tau[ii+1]*cos(thetas[jj]),(a_b*b_b/sqrt((b_b*cos(thetas[jj])**2)+((a_b*sin(thetas[jj]))**2)))*tau[ii+1]*sin(thetas[jj])])
                p3=numpy.array([0,(a_b*b_b/sqrt((b_b*cos(thetas[jj+1])**2)+((a_b*sin(thetas[jj+1]))**2)))*tau[ii+1]*cos(thetas[jj+1]),(a_b*b_b/sqrt((b_b*cos(thetas[jj+1])**2)+((a_b*sin(thetas[jj+1]))**2)))*tau[ii+1]*sin(thetas[jj+1])])
                p4=numpy.array([0,(a_b*b_b/sqrt((b_b*cos(thetas[jj+1])**2)+((a_b*sin(thetas[jj+1]))**2)))*tau[ii]*cos(thetas[jj+1]),(a_b*b_b/sqrt((b_b*cos(thetas[jj+1])**2)+((a_b*sin(thetas[jj+1]))**2)))*tau[ii]*sin(thetas[jj+1])])

                n_cap=numpy.array([1,0,0])
                dA=0.5*(numpy.linalg.norm(numpy.cross(p2-p1,p1-p3))+numpy.linalg.norm(numpy.cross(p3-p1,p1-p4)))
            else:
                p1=numpy.array([0,(a_b*b_b/sqrt((b_b*cos(thetas[jj])**2)+((a_b*sin(thetas[jj]))**2)))*tau[ii]*cos(thetas[jj]),(a_b*b_b/sqrt((b_b*cos(thetas[jj])**2)+((a_b*sin(thetas[jj]))**2)))*tau[ii]*sin(thetas[jj])])
                p2=numpy.array([0,(a_b*b_b/sqrt((b_b*cos(thetas[jj])**2)+((a_b*sin(thetas[jj]))**2)))*tau[ii+1]*cos(thetas[jj]),(a_b*b_b/sqrt((b_b*cos(thetas[jj])**2)+((a_b*sin(thetas[jj]))**2)))*tau[ii+1]*sin(thetas[jj])])
                p3=numpy.array([0,(a_b*b_b/sqrt((b_b*cos(thetas[jj+1])**2)+((a_b*sin(thetas[jj+1]))**2)))*tau[ii+1]*cos(thetas[jj+1]),(a_b*b_b/sqrt((b_b*cos(thetas[jj+1])**2)+((a_b*sin(thetas[jj+1]))**2)))*tau[ii+1]*sin(thetas[jj+1])])

                n_cap=numpy.array([1,0,0])
                dA=0.5*numpy.linalg.norm(numpy.cross(p2-p1,p1-p3))

            Pe=P
            Pedot_Y=Pdot_Y

            PT=Pe
            PTdot_Y=Pedot_Y
            PTdot_c=0

            Force=PT*dA*n_cap
            Forcedot_Y=PTdot_Y*dA*n_cap
            Forcedot_c=PTdot_c*dA*n_cap

            F=F+Force
            Fdot_Y=Fdot_Y+Forcedot_Y
            Fdot_c=Fdot_c+Forcedot_c

    fx=2*F[0]
    fz=2*F[2]

    fxdot_Y=2*Fdot_Y[0]
    fzdot_Y=2*Fdot_Y[2]

    fxdot_V=2*Fdot_V[0]
    fzdot_V=2*Fdot_V[2]

    fxdot_c=2*Fdot_c[0]
    fzdot_c=2*Fdot_c[2]

    aerocoeff = numpy.array([fx,fz,fxdot_Y,fzdot_Y,fxdot_V,fzdot_V,fxdot_c,fzdot_c])
    return aerocoeff

def aero_autodiff_u(V,Y,c,alpha):
    a_b=0.5
    b_b=0.5

    n=11

    v0 = 3000.0
    y0 = 40000.0

    v=V*v0

    Rgas=287.0
    rho0 = 1.225
    HScale = 8500.0
    rho=rho0*exp(-(Y*y0/HScale))

    P=101325*exp(-(Y*y0/HScale))

    T=P/(rho*Rgas)

    gamma=1.4

    Mach=v/sqrt(gamma*Rgas*T)

    epsi=(((gamma-1)*(Mach**2))+2)/((gamma+1)*(Mach**2))

    W=numpy.array([-cos(alpha),0,-sin(alpha)])
    Wdot_alpha = numpy.array([sin(alpha),0,-cos(alpha)])

    F=numpy.array([0,0,0])
    Fdot_alpha=numpy.array([0,0,0])

    x=numpy.linspace(0,c,n)
    thetas=numpy.linspace(pi/2,3*pi/2,n)

    for ii in range(n-1):
        for jj in range(n-1):
            if ii!=n-2:
                X=x[ii]

                Xp=x[ii+1]

                thetass=thetas[jj]
                thetap=thetas[jj+1]

                p1=numpy.array([X,a_b*sqrt(1-(X/c))*cos(thetass),b_b*sqrt(1-(X/c))*sin(thetass)])

                p2=numpy.array([Xp,a_b*sqrt(1-(Xp/c))*cos(thetass),b_b*sqrt(1-(Xp/c))*sin(thetass)])

                p3=numpy.array([Xp,a_b*sqrt(1-(Xp/c))*cos(thetap),b_b*sqrt(1-(Xp/c))*sin(thetap)])

                p4=numpy.array([X,a_b*sqrt(1-(X/c))*cos(thetap),b_b*sqrt(1-(X/c))*sin(thetap)])

                N=numpy.cross(p2-p1,p3-p1)

                n_cap=N/numpy.linalg.norm(N)

                cross1 = numpy.cross(p2-p1,p1-p3)
                normcross1 = numpy.linalg.norm(cross1)

                cross2 = numpy.cross(p3-p1,p1-p4)
                normcross2 = numpy.linalg.norm(cross2)

                dA = 0.5*(normcross1 + normcross2)
            else:
                X=x[ii]

                Xp=x[ii+1]

                thetass=thetas[jj]
                thetap=thetas[jj+1]

                p1=numpy.array([X,a_b*sqrt(1-(X/c))*cos(thetass),b_b*sqrt(1-(X/c))*sin(thetass)])

                p2=numpy.array([Xp,a_b*sqrt(1-(Xp/c))*cos(thetass),b_b*sqrt(1-(Xp/c))*sin(thetass)])

                p3=numpy.array([X,a_b*sqrt(1-(X/c))*cos(thetap),b_b*sqrt(1-(X/c))*sin(thetap)])

                N=numpy.cross(p2-p1,p3-p1)

                normN = numpy.linalg.norm(N)

                n_cap=N/normN

                cross1 = numpy.cross(p2-p1,p1-p3)

                normcross1 = numpy.linalg.norm(cross1)

                dA=0.5*normcross1

            costheta=numpy.dot(W,n_cap)
            costhetadot_alpha = numpy.dot(Wdot_alpha,n_cap)

            if costheta<0:
                Pe=P
                Pedot_alpha = 0
            else:
                Pe=P+(0.5*rho*v*v*(2-epsi)*(costheta**2))
                Pedot_alpha=rho*v*v*(2-epsi)*costheta*costhetadot_alpha

            PT=Pe
            PTdot_alpha=Pedot_alpha

            Force=PT*dA*n_cap
            Forcedot_alpha=PTdot_alpha*dA*n_cap

            F=F+Force
            Fdot_alpha=Fdot_alpha+Forcedot_alpha

    # Base
    n=11
    tau=numpy.linspace(0,1,n);

    for ii in range(n-1):
        for jj in range(n-1):
            if ii!=0:
                p1=numpy.array([0,(a_b*b_b/sqrt((b_b*cos(thetas[jj])**2)+((a_b*sin(thetas[jj]))**2)))*tau[ii]*cos(thetas[jj]),(a_b*b_b/sqrt((b_b*cos(thetas[jj])**2)+((a_b*sin(thetas[jj]))**2)))*tau[ii]*sin(thetas[jj])])
                p2=numpy.array([0,(a_b*b_b/sqrt((b_b*cos(thetas[jj])**2)+((a_b*sin(thetas[jj]))**2)))*tau[ii+1]*cos(thetas[jj]),(a_b*b_b/sqrt((b_b*cos(thetas[jj])**2)+((a_b*sin(thetas[jj]))**2)))*tau[ii+1]*sin(thetas[jj])])
                p3=numpy.array([0,(a_b*b_b/sqrt((b_b*cos(thetas[jj+1])**2)+((a_b*sin(thetas[jj+1]))**2)))*tau[ii+1]*cos(thetas[jj+1]),(a_b*b_b/sqrt((b_b*cos(thetas[jj+1])**2)+((a_b*sin(thetas[jj+1]))**2)))*tau[ii+1]*sin(thetas[jj+1])])
                p4=numpy.array([0,(a_b*b_b/sqrt((b_b*cos(thetas[jj+1])**2)+((a_b*sin(thetas[jj+1]))**2)))*tau[ii]*cos(thetas[jj+1]),(a_b*b_b/sqrt((b_b*cos(thetas[jj+1])**2)+((a_b*sin(thetas[jj+1]))**2)))*tau[ii]*sin(thetas[jj+1])])

                n_cap=numpy.array([1,0,0])
                dA=0.5*(numpy.linalg.norm(numpy.cross(p2-p1,p1-p3))+numpy.linalg.norm(numpy.cross(p3-p1,p1-p4)))
            else:
                p1=numpy.array([0,(a_b*b_b/sqrt((b_b*cos(thetas[jj])**2)+((a_b*sin(thetas[jj]))**2)))*tau[ii]*cos(thetas[jj]),(a_b*b_b/sqrt((b_b*cos(thetas[jj])**2)+((a_b*sin(thetas[jj]))**2)))*tau[ii]*sin(thetas[jj])])
                p2=numpy.array([0,(a_b*b_b/sqrt((b_b*cos(thetas[jj])**2)+((a_b*sin(thetas[jj]))**2)))*tau[ii+1]*cos(thetas[jj]),(a_b*b_b/sqrt((b_b*cos(thetas[jj])**2)+((a_b*sin(thetas[jj]))**2)))*tau[ii+1]*sin(thetas[jj])])
                p3=numpy.array([0,(a_b*b_b/sqrt((b_b*cos(thetas[jj+1])**2)+((a_b*sin(thetas[jj+1]))**2)))*tau[ii+1]*cos(thetas[jj+1]),(a_b*b_b/sqrt((b_b*cos(thetas[jj+1])**2)+((a_b*sin(thetas[jj+1]))**2)))*tau[ii+1]*sin(thetas[jj+1])])

                n_cap=numpy.array([1,0,0])
                dA=0.5*numpy.linalg.norm(numpy.cross(p2-p1,p1-p3))

            Pe=P

            PT=Pe

            Force=PT*dA*n_cap

            F=F+Force

    fx=2*F[0]
    fz=2*F[2]

    fxdot_alpha=2*Fdot_alpha[0]
    fzdot_alpha=2*Fdot_alpha[2]

    aerocoeff = numpy.array([fx,fz,fxdot_alpha,fzdot_alpha])
    return aerocoeff
