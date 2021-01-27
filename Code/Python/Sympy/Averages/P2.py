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
    return a,b,c,(10*qaj**2/3 - 31*qaj*qajp1/3 + 11*qaj*qajp2/3 + 25*qajp1**2/3 - 19*qajp1*qajp2/3 + 4*qajp2**2/3)

def P2jm1tojp1(qajp1,qaj,qajm1,dx):
    a =-qaj/dx**2 + qajm1/(2*dx**2) + qajp1/(2*dx**2)
    b=-qajm1/(2*dx) + qajp1/(2*dx)
    c=13*qaj/12 - qajm1/24 - qajp1/24
    return a,b,c,(13*qaj**2/3 - 13*qaj*qajm1/3 - 13*qaj*qajp1/3 + 4*qajm1**2/3 + 5*qajm1*qajp1/3 + 4*qajp1**2/3)


def P2jm2toj(qajm2,qajm1,qaj,dx):
    a = (qaj - 2*qajm1 + qajm2)/ (2*dx**2)
    b = (qajm2 - 4*qajm1 + 3*qaj)/ (2*dx)
    c =  -qajm2/24 + qajm1/12 + 23*qaj/24
    return a,b,c,(10*qaj**2/3 - 31*qaj*qajm1/3 + 11*qaj*qajm2/3 + 25*qajm1**2/3 - 19*qajm1*qajm2/3 + 4*qajm2**2/3)

def P4jm2tojp2(qajm2,qajm1,qaj,qajp1,qajp2,dx):
    a =  qaj/(4*dx**4) - qajm1/(6*dx**4) + qajm2/(24*dx**4) - qajp1/(6*dx**4) + qajp2/(24*dx**4)
    b = qajm1/(6*dx**3) - qajm2/(12*dx**3) - qajp1/(6*dx**3) + qajp2/(12*dx**3)
    c =-11*qaj/(8*dx**2) + 3*qajm1/(4*dx**2) - qajm2/(16*dx**2) + 3*qajp1/(4*dx**2) - qajp2/(16*dx**2)
    d = -17*qajm1/(24*dx) + 5*qajm2/(48*dx) + 17*qajp1/(24*dx) - 5*qajp2/(48*dx)
    e = 1067*qaj/960 - 29*qajm1/480 + 3*qajm2/640 - 29*qajp1/480 + 3*qajp2/640

    
    SI = 77051*qaj**2/1680 - 24923*qaj*qajm1/420 + 7547*qaj*qajm2/560 - 24923*qaj*qajp1/420 + 7547*qaj*qajp2/560 + 104963*qajm1**2/5040 - 51001*qajm1*qajm2/5040 + 89549*qajm1*qajp1/2520 - 38947*qajm1*qajp2/5040 + 1727*qajm2**2/1260 - 38947*qajm2*qajp1/5040 + 8209*qajm2*qajp2/5040 + 104963*qajp1**2/5040 - 51001*qajp1*qajp2/5040 + 1727*qajp2**2/1260
    
    return a,b,c,d,e,SI

#P2 Poly
P2jtojp2a22,P2jtojp2a21,P2jtojp2a20,P2jtojp2aSI = P2jtojp2(qaj,qajp1,qajp2,dx)
P2jm1tojp1a22,P2jm1tojp1a21,P2jm1tojp1a20,P2jm1tojp1aSI = P2jm1tojp1(qajp1,qaj,qajm1,dx)
P2jm2toja22,P2jm2toja21,P2jm2toja20,P2jm2tojaSI = P2jm2toj(qajm2,qajm1,qaj,dx)


#P4 Poly
P4jm2tojp2a44,P4jm2tojp2a43,P4jm2tojp2a42,P4jm2tojp2a41,P4jm2tojp2a40,P4jm2tojp2aSI =P4jm2tojp2(qajm2,qajm1,qaj,qajp1,qajp2,dx)


#j-2 - a2 coeffs
Cjm2a2 = -1/8
Cjm1p1a2 = 5/4
Cjp2a2 = -1/8
RemCoeffa2 = P4jm2tojp2a42 - (-P2jm2toja22/8 + 5*P2jm1tojp1a22/4 - P2jtojp2a22/8)
RemCoeffa2 = RemCoeffa2.simplify()

#a1  coeffs
Cjm2a1 = 5/24
Cjm1p1a1 = 14/24
Cjp2a1 = 5/24
RemCoeffa1 = P4jm2tojp2a41 - (5*P2jm2toja21/24 + 14*P2jm1tojp1a21/24 + 5*P2jtojp2a21/24)
RemCoeffa1 = RemCoeffa1.simplify()


#a0  coeffs
Cjm2a0 = -9/80
Cjm1p1a0 = 98/80
Cjp2a0 = -9/80
RemCoeffa0 = P4jm2tojp2a40 - (-9*P2jm2toja20/80 + 98*P2jm1tojp1a20/80 -9*P2jtojp2a20/80)
RemCoeffa0 = RemCoeffa0.simplify()

# Coeffa2 = linsolve([P4jm2tojp2a42 - a*P2jtojp2a22 -b*P2jm1tojp1a22 - c*P2jm2toja22], (a,b,c))
# Coeffa1 = linsolve([P4jm2tojp2a41 - a*P2jtojp2a21 -b*P2jm1tojp1a21 - c*P2jm2toja21], (a,b,c))
# Coeffa0 = linsolve([P4jm2tojp2a40 - a*P2jtojp2a20 -b*P2jm1tojp1a20 - c*P2jm2toja20], (a,b,c))