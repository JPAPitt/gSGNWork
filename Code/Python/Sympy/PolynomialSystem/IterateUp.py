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

def GeneratePoly(CL,ub,lb,xmxj,qa):
    r = len(CL)
    Poly  = qa
    
    for i in range(r):
        Term = CL[i]*xmxj**(i+1)
        Corr = integrate( Term,(xmxj,lb,ub))/dx
        Poly = Poly + (Term - Corr)
    Poly = Poly.expand()
    Poly = collect(Poly,xmxj)
    return Poly
    
    
#Symbol definition
qaj,dx = symbols('qaj dx')

#Coeffs
a1,a2,a3 = symbols('a1 a2 a3')
xmxj= symbols('xmxj')

#Poly
# print('0')
# Poly0  = GeneratePoly([],dx/2,-dx/2,xmxj,qaj)
# print('1')
# Poly1  = GeneratePoly([a1],dx/2,-dx/2,xmxj,qaj)

# print('2')
# Poly2jjm1  = GeneratePoly([a1,a2],dx/2,-3*dx/2,xmxj,qaj)


# print('3')
# Poly3  = GeneratePoly([a1,a2,a3],dx/2,-dx/2,xmxj,qaj)



