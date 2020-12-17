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

def GenPoly(k,r,ci,qal,cix,x,xl,dx,xj,xmxj):
    
    PolyC= GeneratePoly(k,r,ci,qal,cix,x,xl,dx)
    PolyC= PolyC.subs(x - xj, xmxj)
    PolyC= PolyC.subs(2*x - 2*xj, 2*xmxj)
    PolyC = PolyC.expand()
    PolyC = collect(PolyC,xmxj)
    
    PolyCoeffs = []
    for i in range(k):
        PolyCoeffs.append(PolyC.coeff(xmxj,i))
        
    # print(PolyC)
        
    B= 0
    for l in range(1,k+1):
        Bterm = dx**(2*l - 1) *integrate( diff(PolyC,(xmxj,l))**2, (xmxj,-dx/2,dx/2))
        B = B +  Bterm
    B = B.expand()
    
    return PolyC,PolyCoeffs,B

def ReturnTerm(Poly,xmxj,lb,up,l):
    return integrate( diff(Poly,(xmxj,l)),(xmxj,lb,up))/dx
    
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

#Quad
P4jm2tojp2,P4jm2tojp2Coefficients, SIjP4jm2tojp2 =  GenPoly(5,2,ci,qal,cix,x,xl,dx,xj,xmxj)

#All quadratics from j-2 to j+2
P3jtojp3,P3jtojp3Coefficients, SIjtojp3 =  GenPoly(4,0,ci,qal,cix,x,xl,dx,xj,xmxj)
P3jm1tojp2,P3jm1tojp2Coefficients, SIjm1tojp2 =  GenPoly(4,1,ci,qal,cix,x,xl,dx,xj,xmxj)
P3jm3toj,P3jm3tojCoefficients, SIjm3toj =  GenPoly(4,3,ci,qal,cix,x,xl,dx,xj,xmxj)

#All quadratics from j-2 to j+2
P2jtojp2,P2jtojp2Coefficients, SIjtojp2 =  GenPoly(3,0,ci,qal,cix,x,xl,dx,xj,xmxj)
P2jm2toj,P2jm2tojCoefficients, SIjm2toj =  GenPoly(3,2,ci,qal,cix,x,xl,dx,xj,xmxj)

P2jm1tojp1,P2jm1tojp1Coefficients, SIjm1tojp1 =  GenPoly(3,1,ci,qal,cix,x,xl,dx,xj,xmxj)


#All linears from j-1 to j+1
P1jtojp1,P1jtojp1Coefficients, SIjtojp1 =  GenPoly(2,0,ci,qal,cix,x,xl,dx,xj,xmxj)
P1jm1toj,P1jm1tojCoefficients, SIjm1toj =  GenPoly(2,1,ci,qal,cix,x,xl,dx,xj,xmxj)

#It works!
RegInterpjtojp3 = P3jtojp3
BasisInterpjtojp1 = P1jtojp1
BasisInterpjtojp2 = P2jtojp2Coefficients[-1]*(xmxj**2 - dx*xmxj - dx**2/12)
BasisInterpjtojp3 = P3jtojp3Coefficients[-1]*(xmxj**3 - 3*dx*xmxj**2 + 7*dx**2/4*xmxj + dx**3/4) 
BasisInterpjtojp3Full = P1jtojp1 +  BasisInterpjtojp2 + BasisInterpjtojp3
DiffCheckA = (BasisInterpjtojp3Full - RegInterpjtojp3).expand()

BasisInterpjtojp1C = P1jtojp1
BasisInterpjtojp2C2 = P2jtojp2Coefficients[-1]*(xmxj**2 - dx*xmxj - dx**2/12)
BasisInterpjtojp3C = P3jm1tojp2Coefficients[-1]*(xmxj**3 - 3*dx*xmxj**2 + 7*dx**2/4*xmxj + dx**3/4) 
BasisInterpjtojp3C2Full = BasisInterpjtojp1C + BasisInterpjtojp2C2 + BasisInterpjtojp3C
DiffCheckB = (BasisInterpjtojp3C2Full - RegInterpjtojp3).expand()
# RegInterpjm3toj = P3jm3toj
# BasisInterpjm3toj = P1jm1toj +  P2jm2tojCoefficients[-1]*(xmxj**2 + dx*xmxj - dx**2/12)  + P3jm3tojCoefficients[-1]*(xmxj**3 + 3*dx*xmxj**2 + 7*dx**2/4*xmxj - dx**3/4) 

# BasisInterpjtojp3C = 

# #Quadratic
# #
# Qjm1tojp1 = P2jm1tojp1
# BVjtojp1 = P1jtojp1 +  P2jm1tojp1Coefficients[-1]*(xmxj**2 - dx*xmxj - dx**2/12)
# BVjm1toj = P1jm1toj +  P2jm1tojp1Coefficients[-1]*(xmxj**2 + dx*xmxj - dx**2/12)