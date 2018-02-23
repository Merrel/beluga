from math import *
import numpy

def odes(tau, X, params, aux):

    xCont = aux[0]
    vCont = aux[1]
    flag = aux[2]

    x = X[0]
    v = X[1]
    lamx = X[2]
    lamv = X[3]
    tf = X[4]

    t = tau*tf

    m = 1

    F = -(lamv/(2*m))

    xdot = v
    vdot = F/m
    lamxdot = 0
    lamvdot = -lamx
    tfdot = 0

    XDot = tf*numpy.array([xdot,vdot,lamxdot,lamvdot,tfdot])
    return XDot

def bcs(Ya, Yb, params, aux):
    xCont = aux[0]
    vCont = aux[1]
    flag = aux[2]

    x0 = Ya[0]
    v0 = Ya[1]
    lamx0 = Ya[2]
    lamv0 = Ya[3]
    tf0 = Ya[4]

    x = Yb[0]
    v = Yb[1]
    lamx = Yb[2]
    lamv = Yb[3]
    tf = Yb[4]

    Xdot = odes(1, Yb, params, aux)

    m = 1
    F = -(lamv/(2*m))

    coState = numpy.array([lamx, lamv])

    Ham = 1 + (F**2) + (numpy.dot(coState,Xdot[:2])/tf)

    if flag == 1:
        zeroVec = numpy.array([x0 - 0,\
        v0 - 0,\
        x - xCont,\
        lamv - 0,\
        Ham - 0])
    else:
        zeroVec = numpy.array([x0 - 0,\
        v0 - 0,\
        x - xCont,\
        v - vCont,\
        Ham - 0])

    return zeroVec
