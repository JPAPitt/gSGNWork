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
   
    
#Symbol definition
qaj,qjmh,qjph,dx = symbols('qaj qjmh qjph dx')
x,xj,xmxj,dx = symbols('x xj xmxj dx')

pa,pb,pc = symbols('pa pb pc')

Poly = pa*(xmxj)**2 + pb*xmxj + pc

PolyAvgjmhjph = integrate(Poly,(xmxj,-dx/2,dx/2))/dx
Polyatjmh = Poly.subs(xmxj,-dx/2)
Polyatjph = Poly.subs(xmxj,dx/2)

PolyThroughPointsAndAverage = linsolve([PolyAvgjmhjph - qaj, Polyatjmh - qjmh, Polyatjph - qjph], (pa,pb,pc))


