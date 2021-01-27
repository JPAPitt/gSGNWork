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

def P1jtojp1(qaj,qajp1,dx):
    return (qajp1 - qaj)/ dx ,qaj,(qaj**2 - 2*qaj*qajp1 + qajp1**2)

def P1jm1toj(qaj,qajm1,dx):
    return (qaj - qajm1)/ dx ,qaj,(qaj**2 - 2*qaj*qajm1 + qajm1**2)


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

#P1 Poly
P1jtojp1a11,P1jtojp1a10,P1jtojp1SI  =  P1jtojp1(qaj,qajp1,dx)
P1jm1toja11,P1jm1toja10,P1jm1tojSI  =  P1jm1toj(qaj,qajm1,dx)

P1jm2tojm1a11,P1jm2tojm1a10,P1jm2tojm1SI  =  P1jm1toj(qajm1,qajm2,dx)
P1jp1tojp2a11,P1jp1tojp2a10,P1jp1tojp2SI  =  P1jtojp1(qajp1,qajp2,dx)

#P2 Polys
P2jm2toja22,P2jm2toja21,P2jm2toja20,P2jm2tojSI = P2jm2toj(qajm2,qajm1,qaj,dx)   
P2jm1tojp1a22,P2jm1tojp1a21,P2jm1tojp1a20,P2jm1tojp1SI  = P2jm1tojp1(qajp1,qaj,qajm1,dx)
P2jtojp2a22,P2jtojp2a21,P2jtojp2a20,P2jtojp2SI =   P2jtojp2(qaj,qajp1,qajp2,dx) 



#Create P2jm1tojp1a21,P2jm1tojp1a20 From P1jm1toja11,P1jm1toja10 and P1jtojp1a11,P1jtojp1a10 


CPjm1ja1 = 1/2
CPjjp1a1 = 1/2
RemCoeffPjm1jp1a1 = P2jm1tojp1a21 - (P1jtojp1a11/2 + P1jm1toja11/2)
RemCoeffPjm1jp1a1 = RemCoeffPjm1jp1a1.simplify()

CPjm1ja0 = 1/2
CPjjp1a0 = 1/2
RemCoeffPjm1jp1a0 = P2jm1tojp1a20 - (P1jtojp1a10/2 + P1jm1toja10/2)
RemCoeffPjm1jp1a0 = RemCoeffPjm1jp1a0.simplify()



# # #Create P2jm2toja21,P2jm2toja20 From P1jm1toja11,P1jm1toja10 and P1jm2tojm1a11,P1jm2tojm1a10 
# CPjm2jm1a1 = -1/2
# CPjm1ja1 = 3/2
# RemCoeffPjm2ja1 = P2jm2toja21 - (-P1jm2tojm1a11/2 + 3*P1jm1toja11/2)
# RemCoeffPjm2ja1 = RemCoeffPjm2ja1.simplify()



# # #Create P2jtojp2a21,P2jtojp2a20 From P1jp1tojp2a11,P1jp1tojp2a10 and P1jtojp1a11,P1jtojp1a10
# CPjp1jp2a1 = -1/2
# CPjjp1a1 = 3/2
# RemCoeffPjjp2a1 =P2jtojp2a21 - (-P1jp1tojp2a11/2 + 3*P1jtojp1a11/2)
# RemCoeffPjjp2a1 =RemCoeffPjjp2a1.simplify()




#Create P2jm2toja21,P2jm2toja20 From P1jm1toja11,P1jm1toja10 and P1jm2tojm1a11,P1jm2tojm1a10 with appropriate power of dx
# CPjm2jm1a0 = 1/12
# CPjm1ja0 = 23/24
# RemCoeffPjm2ja0 = P2jm2toja20 - (P1jm2tojm1a10/12 + 11*P1jm1toja10/12)
# RemCoeffPjm2ja0 = RemCoeffPjm2ja0.simplify()

# a0qajm2dx = qajm2




# #j-2 - a2 coeffs
# Cjm2a2 = -1/8
# Cjm1p1a2 = 5/4
# Cjp2a2 = -1/8
# RemCoeffa2 = P4jm2tojp2a42 - (-P2jm2toja22/8 + 5*P2jm1tojp1a22/4 - P2jtojp2a22/8)
# RemCoeffa2 = RemCoeffa2.simplify()

# #a1  coeffs
# Cjm2a1 = 5/24
# Cjm1p1a1 = 14/24
# Cjp2a1 = 5/24
# RemCoeffa1 = P4jm2tojp2a41 - (5*P2jm2toja21/24 + 14*P2jm1tojp1a21/24 + 5*P2jtojp2a21/24)
# RemCoeffa1 = RemCoeffa1.simplify()


# #a0  coeffs
# Cjm2a0 = 9/80
# Cjm1p1a0 = 62/80
# Cjp2a0 = 9/80
# RemCoeffa0 = P4jm2tojp2a40 - (9*P2jm2toja21/80 + 62*P2jm1tojp1a21/80 + 9*P2jtojp2a21/80)
# RemCoeffa0 = RemCoeffa1.simplify()

# # Coeffa2 = linsolve([P4jm2tojp2a42 - a*P2jtojp2a22 -b*P2jm1tojp1a22 - c*P2jm2toja22], (a,b,c))
# # Coeffa1 = linsolve([P4jm2tojp2a41 - a*P2jtojp2a21 -b*P2jm1tojp1a21 - c*P2jm2toja21], (a,b,c))
# # Coeffa0 = linsolve([P4jm2tojp2a40 - a*P2jtojp2a20 -b*P2jm1tojp1a20 - c*P2jm2toja20], (a,b,c))