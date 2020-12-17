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
a0,a1,a2,a3 = symbols('a0 a1 a2 a3')
b2 = symbols('b2')
xmxj= symbols('xmxj')

Poly0 = a0
Poly1 = a1*xmxj + a0
Poly2 = a2*xmxj**2 + a1*xmxj + a0

Basis0 = a0
Basis1 = a1*xmxj - integrate(a1*xmxj,(xmxj,-dx/2,dx/2))/dx
Basis2jm1j = a2*xmxj**2 - integrate(a2*xmxj,(xmxj,-3*dx/2,-dx/2))/(dx) *xmxj  - integrate(a2*xmxj**2,(xmxj,-dx/2,dx/2))/dx
Basis2jjp1 = a2*xmxj**2 - integrate(a2*xmxj,(xmxj,dx/2,3*dx/2))/(dx) *xmxj - integrate(a2*xmxj**2,(xmxj,-dx/2,dx/2))/dx

# Basis3jm2jm1j= a3*xmxj**3  - integrate(a3*xmxj,(xmxj,-3*dx/2,-dx/2))/(dx) *xmxj   - integrate(a3*xmxj**3,(xmxj,-dx/2,dx/2))/dx
Basis3jm2jm1jCA = a3*xmxj**3   - integrate(a3*xmxj**3,(xmxj,-dx/2,dx/2))/dx
Basis3jm2jm1jZjm1 = integrate(a3*xmxj**2,(xmxj,-3*dx/2,-dx/2))/(dx) *xmxj
Basis3jm2jm1jZjm2 = integrate(a3*xmxj,(xmxj,-5*dx/2,-3*dx/2))/(dx) *xmxj**2
#- integrate(a3*xmxj**2,(xmxj,-3*dx/2,-dx/2))/(dx) *xmxj

# Basis3jm2tojIj = integrate(Basis3jm2jm1j,(xmxj,-dx/2,dx/2))/dx
# Basis3jm2tojIjm1 = integrate(Basis3jm2jm1j,(xmxj,-3*dx/2,-dx/2))/dx

Basis3jm2jm1jCAIjm1 = integrate(Basis3jm2jm1jCA,(xmxj,-3*dx/2,-dx/2))/dx
Basis3jm2jm1jZjm1Ijm1 = integrate(Basis3jm2jm1jZjm1,(xmxj,-3*dx/2,-dx/2))/dx
Basis3jm2jm1jZjm1Ijm2 = integrate(Basis3jm2jm1jZjm1,(xmxj,-5*dx/2,-3*dx/2))/dx

# Basis2jjp1CA = a2*xmxj**2 - integrate(a2*xmxj**2,(xmxj,-dx/2,dx/2))/dx
# Basis2jjp1Lin =integrate(a2*xmxj,(xmxj,dx/2,3*dx/2))/(dx) *xmxj 

# Basis1Ij = integrate(Basis1,(xmxj,-dx/2,dx/2))/dx
# Basis1Ijp1 = integrate(Basis1,(xmxj,dx/2,3*dx/2))/dx
# Basis2jjp1Ij = integrate(Basis2jjp1,(xmxj,-dx/2,dx/2))/dx
# Basis2jjp1Ijp1 = integrate(Basis2jjp1,(xmxj,dx/2,3*dx/2))/dx

# Basis2jm1jIj = integrate(Basis2jm1j,(xmxj,-dx/2,dx/2))/dx
# Basis2jm1jm1 = integrate(Basis2jm1j,(xmxj,-3*dx/2,-dx/2))/dx
# Basis2jjp1CAIj = integrate(Basis2jjp1CA ,(xmxj,-dx/2,dx/2))/dx
# Basis2jjp1CAIjp1 = integrate(Basis2jjp1CA,(xmxj,dx/2,3*dx/2))/dx

# Basis2jjp1LinIj = integrate(Basis2jjp1Lin ,(xmxj,-dx/2,dx/2))/dx
# Basis2jjp1LinIjp1 = integrate(Basis2jjp1Lin,(xmxj,dx/2,3*dx/2))/dx

# Basis2Ijm1 = integrate(Basis2jm1j,(xmxj,-3*dx/2,-dx/2))/dx
# Poly0I = integrate(Poly0,(xmxj,-dx/2,dx/2))/dx

#Poly
# print('0')
# Poly0  = GeneratePoly([],dx/2,-dx/2,xmxj,qaj)
# print('1')
# Poly1  = GeneratePoly([a1],dx/2,-dx/2,xmxj,qaj)

# print('2')
# Poly2jjm1  = GeneratePoly([a1,a2],dx/2,-3*dx/2,xmxj,qaj)


# print('3')
# Poly3  = GeneratePoly([a1,a2,a3],dx/2,-dx/2,xmxj,qaj)



