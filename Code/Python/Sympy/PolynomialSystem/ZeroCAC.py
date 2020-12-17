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

a0,a1,a2,a3,a4,a5,a6 = symbols('a0 a1 a2 a3 a4 a5 a6')
c0,c1,c2,c3,c4,c5,c6 = symbols('c0 c1 c2 c3 c4 c5 c6')
xmxj,dx= symbols('xmxj dx')

#Linear Term
P1j = xmxj - c0
P1jInt = integrate(P1j,(xmxj,-dx/2,dx/2))
P1jCorrj= linsolve([P1jInt], c0)

#Quadratic Term
P2j = xmxj**2 + c1*xmxj + c0
P2jInt = integrate(P2j,(xmxj,-dx/2,dx/2))
P2jm1Int = integrate(P2j,(xmxj,-3*dx/2,-dx/2))
P2jp1Int = integrate(P2j,(xmxj,dx/2,3*dx/2))

P2jCorrjm1toj = linsolve([P2jInt,P2jm1Int], (c1,c0))
P2jCorrjtojp1 = linsolve([P2jInt,P2jp1Int], (c1,c0))


# #Cubic Term
P3j = xmxj**3 + c2*xmxj**2  + c1*xmxj + c0
P3jInt = integrate(P3j,(xmxj,-dx/2,dx/2))
P3jm1Int = integrate(P3j,(xmxj,-3*dx/2,-dx/2))
P3jm2Int = integrate(P3j,(xmxj,-5*dx/2,-3*dx/2))
P3jp1Int = integrate(P3j,(xmxj,dx/2,3*dx/2))
P3jp2Int = integrate(P3j,(xmxj,3*dx/2,5*dx/2))

P3jCorrjm2toj = linsolve([P3jInt,P3jm1Int,P3jm2Int], (c2,c1,c0))
P3jCorrjtojp2 = linsolve([P3jInt,P3jp1Int,P3jp2Int], (c2,c1,c0))
P3jCorrjm1tojp1 = linsolve([P3jm1Int,P3jInt,P3jp1Int], (c2,c1,c0))

    
# Quartics
P4j = xmxj**4 - c3*xmxj**3 - c2*xmxj**2  - c1*xmxj - c0
P4jInt = integrate(P4j,(xmxj,-dx/2,dx/2))
P4jm1Int = integrate(P4j,(xmxj,-3*dx/2,-dx/2))
P4jm2Int = integrate(P4j,(xmxj,-5*dx/2,-3*dx/2))
P4jm3Int = integrate(P4j,(xmxj,-7*dx/2,-5*dx/2))

P4jp1Int = integrate(P4j,(xmxj,dx/2,3*dx/2))
P4jp2Int = integrate(P4j,(xmxj,3*dx/2,5*dx/2))
P4jp3Int = integrate(P4j,(xmxj,5*dx/2,7*dx/2))

P4jCorrjm3toj = linsolve([P4jInt,P4jm1Int,P4jm2Int,P4jm3Int], (c3,c2,c1,c0))
P4jCorrjtojp3 = linsolve([P4jInt,P4jp1Int,P4jp2Int,P4jp3Int], (c3,c2,c1,c0))


