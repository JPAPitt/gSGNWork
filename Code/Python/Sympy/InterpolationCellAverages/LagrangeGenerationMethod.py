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


##### Functions #####

# Gen sum does the inner sum in equation (2.19)
#The definitions of k,m and r follow from Shu's notes
# k corresponds to the desired polynomial for k = 1 we get constants, for k=2 we get linears,etc
# m corresponds to what term we are on 
# r is the shift when r = [0,k-1] of the stencil, which can be written assuming we are reconstrucitng over the jth cell  x_{j-r}, x_{j-r+1}, ... x_{j + k-1 - r} 
# ci is used to give the location of the x_{j-1/2} in the list of locations xl (cell edges)
def GenSum(k,m,r,ci,x,xl,dx):
    Sum = 0
    for l in range(k+1):
        if l != m:
            Prod = 1
            for q in range(k+1):
                if q!=l and q!=m:
                    Prod = Prod*(x - xl[ci -r + q])
            Sum = Sum + Prod
    Norm = 1
    for l in range(k+1):
        if l!=m:
             Norm =  Norm*(m-l)*dx
            
    return Sum/Norm
               
# GeneratePoly uses GenSum does to calculate the interpolating polynomial in (2.19)
#qal is the list of cell averages with average over jth cell at index cix 
def GeneratePoly(k,r,ci,qal,cix,x,xl,dx):
    Poly = 0
    for m in range(k+1):
        for j in range(m):
            Poly = Poly +qal[ci+j-r]*dx*GenSum(k,m,r,cix,x,xl,dx)

    return Poly

#Symbol definition
qaj,qajm1,qajm2,qajm3,qajp1,qajp2,qajp3,dx = symbols('qaj qajm1 qajm2 qajm3 qajp1 qajp2  qajp3 dx')
x,xj,xmxj,dx = symbols('x xj xmxj dx')

#These lists contain the relevant information to calculate the polynomial that agrees with the cell averages in over cells {j-r}, {j-r+1}, ... {j + k-1 - r} 
#These lists must contain all the relevant cells information from left to right

#qal - list of the cell averages in j-3,..j+3 (all required to fit every possible cubic) note that the jth cell is in the qal[3]
qal = [qajm3,qajm2,qajm1,qaj,qajp1,qajp2,qajp3]

#xl - list of left boundary of cells j-3,...,j+3 and the right boundary of j+3
xl = [xj - 7*dx/2,xj - 5*dx/2,xj - 3*dx/2,xj - dx/2,xj+dx/2,xj+3*dx/2,xj+5*dx/2,xj+7*dx/2]



#k is desired order of polynomial
#ci is location of cell average of cell j in qal  
#cix is location of left boundary of cell j in xl             
k = 4
ci = 3
cix = 3

#Cubic from j to j+3
P3jtojp3 = GeneratePoly(k,0,ci,qal,cix,x,xl,dx)
P3jtojp3= P3jtojp3.subs(x - xj, xmxj)
P3jtojp3 = P3jtojp3.expand()
P3jtojp3Poly = collect(P3jtojp3,xmxj)
P3jtojp3Coefficients = [P3jtojp3Poly.coeff(xmxj,3),P3jtojp3Poly.coeff(xmxj,2),P3jtojp3Poly.coeff(xmxj,1),P3jtojp3Poly.coeff(xmxj,0)]

#6th order polynomial from j-3 to j+3
P6jm3tojp3 = GeneratePoly(7,3,ci,qal,cix,x,xl,dx)
P6jm3tojp3= P6jm3tojp3.subs(x - xj, xmxj)
P6jm3tojp3 = P6jm3tojp3.expand()
P6jm3tojp3Poly = collect(P6jm3tojp3,xmxj)
P6jm3tojp3Coefficients = [P6jm3tojp3Poly.coeff(xmxj,6),P6jm3tojp3Poly.coeff(xmxj,5),P6jm3tojp3Poly.coeff(xmxj,4),P6jm3tojp3Poly.coeff(xmxj,3),P6jm3tojp3Poly.coeff(xmxj,2),P6jm3tojp3Poly.coeff(xmxj,1),P6jm3tojp3Poly.coeff(xmxj,0)]

#4th order polynomial from j-2 to j+2
P4jm2tojp2 = GeneratePoly(5,2,ci,qal,cix,x,xl,dx)
P4jm2tojp2= P4jm2tojp2.subs(x - xj, xmxj)
P4jm2tojp2 = P4jm2tojp2.expand()
P4jm2tojp2Poly = collect(P4jm2tojp2,xmxj)
P4jm2tojp2Coefficients = [P4jm2tojp2Poly.coeff(xmxj,4),P4jm2tojp2Poly.coeff(xmxj,3),P4jm2tojp2Poly.coeff(xmxj,2),P4jm2tojp2Poly.coeff(xmxj,1),P4jm2tojp2Poly.coeff(xmxj,0)]

#2nd order polynomial from j-1 to j+1
P3jm1tojp1 = GeneratePoly(3,1,ci,qal,cix,x,xl,dx)
P3jm1tojp1= P3jm1tojp1.subs(x - xj, xmxj)
P3jm1tojp1 = P3jm1tojp1.expand()
P3jm1tojp1Poly = collect(P3jm1tojp1,xmxj)
P3jm1tojp1Coefficients = [P3jm1tojp1Poly.coeff(xmxj,2),P3jm1tojp1Poly.coeff(xmxj,1),P3jm1tojp1Poly.coeff(xmxj,0)]


