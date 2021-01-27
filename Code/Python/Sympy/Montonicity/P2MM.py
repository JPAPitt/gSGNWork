"""
This is a script that helps generate the polynomial of the form

P_i(x) = a_n (x - x_i)^n + ... + a_0

That agrees with the cell averages in the ith cell and its chosen neighbours.
For a description see Shu's work explaining WENO and method to generate these. (https://www3.nd.edu/~zxu2/acms60790S13/Shu-WENO-notes.pdf)

Essetially the idea is if we want to reconstruct a function v(x) using cell average values \bar{v}(x), then we can rewrite this problem
as a more standard interpolation prolbem at points by thinking about the function

V(x) = \int_{-\infty}^{x} v(s) ds

So we are intergrating a suitable antiderivative of v, that will when interpolated give us the polynomial that agrees with the cell average values

"""
##### Imports #####
from sympy import *
from sympy.solvers.solveset import linsolve

qaj,qajm1,qajm2,qajp1,qajp2 = symbols('\\bar{q}_j \\bar{q}_{j-1} \\bar{q}_{j-2} \\bar{q}_{j+1} \\bar{q}_{j+2}')
a,b,c = symbols('a b c')
xmxj,dx= symbols('xmxj dx')

def P2jtojp2(qaj,qajp1,qajp2,dx):
    a = (qajp2 - 2*qajp1 + qaj)/ (2*dx**2)
    b = (-qajp2 + 4*qajp1 - 3*qaj)/ (2*dx)
    c =  -qajp2/24 + qajp1/12 + 23*qaj/24
    return a,b,c

def P2jm1tojp1(qajp1,qaj,qajm1,dx):
    a =-qaj/dx**2 + qajm1/(2*dx**2) + qajp1/(2*dx**2)
    b=-qajm1/(2*dx) + qajp1/(2*dx)
    c=13*qaj/12 - qajm1/24 - qajp1/24
    return a,b,c

def P2jm2toj(qajm2,qajm1,qaj,dx):
    a = (qaj - 2*qajm1 + qajm2)/ (2*dx**2)
    b = (qajm2 - 4*qajm1 + 3*qaj)/ (2*dx)
    c =  -qajm2/24 + qajm1/12 + 23*qaj/24
    return a,b,c

def P4jm2tojp2(qajm2,qajm1,qaj,qajp1,qajp2,dx):
    a =  qaj/(4*dx**4) - qajm1/(6*dx**4) + qajm2/(24*dx**4) - qajp1/(6*dx**4) + qajp2/(24*dx**4)
    b = qajm1/(6*dx**3) - qajm2/(12*dx**3) - qajp1/(6*dx**3) + qajp2/(12*dx**3)
    c =-11*qaj/(8*dx**2) + 3*qajm1/(4*dx**2) - qajm2/(16*dx**2) + 3*qajp1/(4*dx**2) - qajp2/(16*dx**2)
    d = -17*qajm1/(24*dx) + 5*qajm2/(48*dx) + 17*qajp1/(24*dx) - 5*qajp2/(48*dx)
    e = 1067*qaj/960 - 29*qajm1/480 + 3*qajm2/640 - 29*qajp1/480 + 3*qajp2/640

    return a,b,c,d,e

#P2 Poly
P2jtojp2a22,P2jtojp2a21,P2jtojp2a20 = P2jtojp2(qaj,qajp1,qajp2,dx)
P2jm1tojp1a22,P2jm1tojp1a21,P2jm1tojp1a20 = P2jm1tojp1(qajp1,qaj,qajm1,dx)
P2jm2toja22,P2jm2toja21,P2jm2toja20 = P2jm2toj(qajm2,qajm1,qaj,dx)


#P4 Poly
P4jm2tojp2a44,P4jm2tojp2a43,P4jm2tojp2a42,P4jm2tojp2a41,P4jm2tojp2a40 =P4jm2tojp2(qajm2,qajm1,qaj,qajp1,qajp2,dx)


Pjtojp2P = P2jtojp2a22*xmxj**2 + P2jtojp2a21*xmxj + P2jtojp2a20 
Pjm1tojp1P = P2jm1tojp1a22*xmxj**2 +  P2jm1tojp1a21*xmxj + P2jm1tojp1a20
Pjm2tojP = P2jm2toja22*xmxj**2 +  P2jm2toja21*xmxj + P2jm2toja20

DPjtojp2P = diff(Pjtojp2P,xmxj)
DPjm1tojp1P = diff(Pjm1tojp1P,xmxj)
DPjm2tojP = diff(Pjm2tojP,xmxj)

Pjtojp2Patjmh = Pjtojp2P.subs(xmxj,-dx/2)
Pjtojp2Patjph = Pjtojp2P.subs(xmxj,dx/2)

Pjm1tojp1Patjmh = Pjm1tojp1P.subs(xmxj,-dx/2)
Pjm1tojp1Patjmh = Pjm1tojp1Patjmh.simplify()
Pjm1tojp1Patjph = Pjm1tojp1P.subs(xmxj,dx/2)
Pjm1tojp1Patjph = Pjm1tojp1Patjph.simplify()

Pjm2tojPatjmh = Pjm2tojP.subs(xmxj,-dx/2)
Pjm2tojPatjph = Pjm2tojP.subs(xmxj,dx/2)

