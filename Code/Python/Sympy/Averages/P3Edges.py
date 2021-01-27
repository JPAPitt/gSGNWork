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

qaj,qajm1,qajm2,qajm3,qajp1,qajp2,qajp3 = symbols('qaj qajm1 qajm2 qajm3 qajp1 qajp2 qajp3')
a,b,c = symbols('a b c')
xmxj,dx= symbols('xmxj dx')


def P3jtojp3(qaj,qajp1,qajp2,qajp3,dx):
    a=(3*qajp1 - 3*qajp2 + qajp3 - qaj)/(6*dx**3)
    b=(-5*qajp1 + 4*qajp2 - qajp3 + 2*qaj)/(2*dx**2)
    c=(69*qajp1 - 33*qajp2 + 7*qajp3 - 43*qaj)/(24*dx)
    d=5*qajp1/24 - qajp2/6 + qajp3/24 + 11*qaj/12
    
    SmoothnessIndicator = 11003*qajp1**2/240 - 8623*qajp1*qajp2/120 \
    + 2321*qajp1*qajp3/120 - 1567*qajp1*qaj/40 + 7043*qajp2**2/240 \
    - 647*qajp2*qajp3/40 + 3521*qajp2*qaj/120 + 547*qajp3**2/240 \
    - 309*qajp3*qaj/40 + 2107*qaj**2/240        

    return a,b,c,d,SmoothnessIndicator


def P3jm1tojp2(qaj,qajm1,qajp1,qajp2,dx):
    a=(-3*qajp1 + qajp2 - qajm1 + 3*qaj)/(6*dx**3)
    b=(qajp1 + qajm1 - 2*qaj)/(2*dx**2)
    c=(27*qajp1 - 5*qajp2 - 7*qajm1 - 15*qaj)/(24*dx)
    d=-qajp1/24 - qajm1/24 + 13*qaj/12
    
    SmoothnessIndicator = 2843*qajp1**2/240 - 821*qajp1*qajp2/120 \
    + 961*qajp1*qajm1/120 - 2983*qajp1*qaj/120 + 89*qajp2**2/80 \
    - 247*qajp2*qajm1/120 + 267*qajp2*qaj/40 + 547*qajm1**2/240 \
    - 1261*qajm1*qaj/120 + 3443*qaj**2/240

    return a,b,c,d,SmoothnessIndicator

def P3jm2tojp1(qaj,qajm1,qajm2,qajp1,dx):
    a=(qajp1 + 3*qajm1 - qajm2 - 3*qaj)/(6*dx**3)
    b=(qajp1 + qajm1 - 2*qaj)/(2*dx**2)
    c=(7*qajp1- 27*qajm1 + 5*qajm2 + 15*qaj)/(24*dx)
    d=-qajp1/24 - qajm1/24 + 13*qaj/12
    
    SmoothnessIndicator = 547*qajp1**2/240 + 961*qajp1*qajm1/120 \
    - 247*qajp1*qajm2/120 - 1261*qajp1*qaj/120 + 2843*qajm1**2/240 \
    - 821*qajm1*qajm2/120 - 2983*qajm1*qaj/120 + 89*qajm2**2/80 \
    + 267*qajm2*qaj/40 + 3443*qaj**2/240

    return a,b,c,d,SmoothnessIndicator


def P3jm3toj(qaj,qajm1,qajm2,qajm3,dx):
    a=(-3*qajm1 + 3*qajm2 - qajm3 + qaj)/(6*dx**3)
    b=(-5*qajm1 + 4*qajm2 - qajm3 + 2*qaj)/(2*dx**2)
    c=(-69*qajm1 + 33*qajm2 - 7*qajm3 + 43*qaj)/(24*dx)
    d=5*qajm1/24 - qajm2/6 + qajm3/24 + 11*qaj/12
    
    SmoothnessIndicator = 11003*qajm1**2/240 - 8623*qajm1*qajm2/120 \
    + 2321*qajm1*qajm3/120 - 1567*qajm1*qaj/40 + 7043*qajm2**2/240 \
    - 647*qajm2*qajm3/40 + 3521*qajm2*qaj/120 + 547*qajm3**2/240 \
    - 309*qajm3*qaj/40 + 2107*qaj**2/240

    return a,b,c,d,SmoothnessIndicator



def P6jm3tojp3(qajm3,qajm2,qajm1,qaj,qajp1,qajp2,qajp3,dx):
    a = -qaj/(36*dx**6) + qajm1/(48*dx**6) - qajm2/(120*dx**6) + qajm3/(720*dx**6) + qajp1/(48*dx**6) - qajp2/(120*dx**6) + qajp3/(720*dx**6)
    b = -qajm1/(48*dx**5) + qajm2/(60*dx**5) - qajm3/(240*dx**5) + qajp1/(48*dx**5) - qajp2/(60*dx**5) + qajp3/(240*dx**5)
    c = 61*qaj/(144*dx**4) - 19*qajm1/(64*dx**4) + 3*qajm2/(32*dx**4) - 5*qajm3/(576*dx**4) - 19*qajp1/(64*dx**4) + 3*qajp2/(32*dx**4) - 5*qajp3/(576*dx**4)
    d =83*qajm1/(288*dx**3) - 13*qajm2/(72*dx**3) + 7*qajm3/(288*dx**3) - 83*qajp1/(288*dx**3) + 13*qajp2/(72*dx**3) - 7*qajp3/(288*dx**3)
    e =-301*qaj/(192*dx**2) + 229*qajm1/(256*dx**2) - 77*qajm2/(640*dx**2) + 37*qajm3/(3840*dx**2) + 229*qajp1/(256*dx**2) - 77*qajp2/(640*dx**2) + 37*qajp3/(3840*dx**2)
    f = -1891*qajm1/(2304*dx) + 559*qajm2/(2880*dx) - 259*qajm3/(11520*dx) + 1891*qajp1/(2304*dx) - 559*qajp2/(2880*dx) + 259*qajp3/(11520*dx)
    g = 30251*qaj/26880 - 7621*qajm1/107520 + 159*qajm2/17920 - 5*qajm3/7168 - 7621*qajp1/107520 + 159*qajp2/17920 - 5*qajp3/7168

    
    SI = 1607739169*qaj**2/2993760 - 790531177*qaj*qajm1/997920 + 1506944981*qaj*qajm2/4989600 - 701563133*qaj*qajm3/14968800 - 790531177*qaj*qajp1/997920 + 1506944981*qaj*qajp2/4989600 - 701563133*qaj*qajp3/14968800 + 16790707*qajm1**2/55440 - 176498513*qajm1*qajm2/739200 + 761142961*qajm1*qajm3/19958400 + 250523543*qajm1*qajp1/443520 - 464678369*qajm1*qajp2/2217600 + 158544319*qajm1*qajp3/4989600 + 108444169*qajm2**2/2217600 - 323333323*qajm2*qajm3/19958400 - 464678369*qajm2*qajp1/2217600 + 84263749*qajm2*qajp2/1108800 - 225623953*qajm2*qajp3/19958400 + 2627203*qajm3**2/1871100 + 158544319*qajm3*qajp1/4989600 - 225623953*qajm3*qajp2/19958400 + 99022657*qajm3*qajp3/59875200 + 16790707*qajp1**2/55440 - 176498513*qajp1*qajp2/739200 + 761142961*qajp1*qajp3/19958400 + 108444169*qajp2**2/2217600 - 323333323*qajp2*qajp3/19958400 + 2627203*qajp3**2/1871100    
    return a,b,c,d,e,f,g,SI


def GetTerm(Exp1,ZeroQList,GetQ):
    Exp1Term = Exp1
    n = len(ZeroQList)
    for i in range(n):
        Exp1Term = Exp1Term.subs(ZeroQList[i],0)
    Exp1Term = Exp1Term.subs(GetQ,1)
    
    return Exp1Term

def CoeffRemTerm(Exp1,Exp2,ZeroQList,GetQ):
    Exp1Term = GetTerm(Exp1,ZeroQList,GetQ)
    Exp2Term = GetTerm(Exp2,ZeroQList,GetQ)
    
    Coeff = Exp1Term / Exp2Term
    
    Rem = Exp1 - Coeff*Exp2
    
    return Coeff,Rem.simplify()

def RepeatCancellation(StartExp,ExpList,AllQList):
    n = len(ExpList)
    CoeffList = [0]*n
    RemExp = StartExp
    for i in range(n-1):
        
        CoeffList[i], RemExp = CoeffRemTerm(RemExp,ExpList[i],AllQList[:i] + AllQList[i+1:],AllQList[i])
    
    CoeffList[n-1] = (RemExp / ExpList[n-1]).simplify()
    
    Check = StartExp
    for i in range(n):
        Check = Check - CoeffList[i]*ExpList[i]
    Check = Check.simplify()
    
    return CoeffList,Check
        

#P3 Poly
P3jm3toja33,P3jm3toja32,P3jm3toja31,P3jm3toja30,P3jm3tojSI = P3jm3toj(qaj,qajm1,qajm2,qajm3,dx)
P3jm2tojp1a33,P3jm2tojp1a32,P3jm2tojp1a31,P3jm2tojp1a30,P3jm2tojp1SI = P3jm2tojp1(qaj,qajm1,qajm2,qajp1,dx)
P3jm1tojp2a33,P3jm1tojp2a32,P3jm1tojp2a31,P3jm1tojp2a30,P3jm1tojp2SI = P3jm1tojp2(qaj,qajm1,qajp1,qajp2,dx)
P3jtojp3a33,P3jtojp3a32,P3jtojp3a31,P3jtojp3a30,P3jtojp3SI = P3jtojp3(qaj,qajp1,qajp2,qajp3,dx)


#P3 Edges
P3jm3tojmh =  P3jm3toja33*(-dx/2)**3  + P3jm3toja32*(-dx/2)**2 + P3jm3toja31*(-dx/2) + P3jm3toja30
P3jm3tojph =  P3jm3toja33*(dx/2)**3  + P3jm3toja32*(dx/2)**2 + P3jm3toja31*(dx/2) + P3jm3toja30


P3jm2tojp1mh = P3jm2tojp1a33*(-dx/2)**3  + P3jm2tojp1a32*(-dx/2)**2  +P3jm2tojp1a31*(-dx/2) +P3jm2tojp1a30
P3jm2tojp1ph = P3jm2tojp1a33*(dx/2)**3  + P3jm2tojp1a32*(dx/2)**2  +P3jm2tojp1a31*(dx/2) +P3jm2tojp1a30

P3jm1tojp2mh = P3jm1tojp2a33*(-dx/2)**3  +P3jm1tojp2a32*(-dx/2)**2  +P3jm1tojp2a31*(-dx/2) +P3jm1tojp2a30
P3jm1tojp2ph = P3jm1tojp2a33*(dx/2)**3  + P3jm1tojp2a32*(dx/2)**2  +P3jm1tojp2a31*(dx/2) +P3jm1tojp2a30

P3jtojp3mh = P3jtojp3a33*(-dx/2)**3  +P3jtojp3a32*(-dx/2)**2  +P3jtojp3a31*(-dx/2) +P3jtojp3a30
P3jtojp3ph = P3jtojp3a33*(dx/2)**3  + P3jtojp3a32*(dx/2)**2  +P3jtojp3a31*(dx/2) +P3jtojp3a30