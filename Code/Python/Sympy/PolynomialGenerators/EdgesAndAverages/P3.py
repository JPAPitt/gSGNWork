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
qaj,qajm1,qajp1,qjmh,qjph,dx = symbols('qaj qajm1 qajp1 qjmh qjph dx')
x,xj,xmxj,dx = symbols('x xj xmxj dx')

pa,pb,pc ,pd= symbols('pa pb pc pd')

Poly = pa*(xmxj)**3 + pb*xmxj**2 + pc*xmxj + pd

PolyAvgjm3hjmh = integrate(Poly,(xmxj,-3*dx/2,-dx/2))/dx
PolyAvgjmhjph = integrate(Poly,(xmxj,-dx/2,dx/2))/dx
PolyAvgjphjp3h = integrate(Poly,(xmxj,dx/2,3*dx/2))/dx
Polyatjmh = Poly.subs(xmxj,-dx/2)
Polyatjph = Poly.subs(xmxj,dx/2)

P3jjm1 = linsolve([PolyAvgjmhjph - qaj, PolyAvgjm3hjmh - qajm1 , Polyatjmh - qjmh, Polyatjph - qjph], (pa,pb,pc,pd))

P3jjp1 = linsolve([PolyAvgjmhjph - qaj, PolyAvgjphjp3h - qajp1, Polyatjmh - qjmh, Polyatjph - qjph], (pa,pb,pc,pd))


