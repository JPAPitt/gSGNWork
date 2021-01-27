"""
This is a python script that tests reconstructions - in particular this one chooses the reconstruction using a minmod comparison of the coefficients of the possible
interpolating cubics. 

This reconstruction - maintains cell average value in cell (due to higher degree terms - their cell average being added) [Easy to prove]
                    - obtains the sufficient order due to use of the required polynomials[Easy to prove]
                    -obtains TVD property since all coefficients will be zero when function in no smooth [Conjecture - not clear that minmod ensures TVD , should be straightforward to show by induction, since the limiting is done recursivley]
                    

My simple test case is this : Piecewise polynomials with sufficient degree (to introduce some error in our apprpoximation - in this case at least 7th) with a discontinuity at x = 0
this allows one to see how discontinuities in higher derivatives are handled by the method, as well as standard jump discontinuities.

"""

##### Imports #####
from numpy import *
from numpy.linalg import solve,norm
from matplotlib.pyplot import plot,loglog,close,legend,title,ylabel,xlabel


##### Interpolation Polynomial Functions #####
"""
This group of functions takes the values of the cell average quantity in the relevant cells (included in abbreviated form in the name)
and returns the coefficient of the cubic that agrees with all the cell averages, as a polynomial of (x - x_j) where j is the cell we are reconstructing.

I have also included the smoothness indicators, used to weight the different coefficients in WENO
"""
#constant
def P0j(q,j):
    return q[j]


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

def P4jm2tojp2(qajm2,qajm1,qaj,qajp1,qajp2,dx):
    a =  qaj/(4*dx**4) - qajm1/(6*dx**4) + qajm2/(24*dx**4) - qajp1/(6*dx**4) + qajp2/(24*dx**4)
    b = qajm1/(6*dx**3) - qajm2/(12*dx**3) - qajp1/(6*dx**3) + qajp2/(12*dx**3)
    c =-11*qaj/(8*dx**2) + 3*qajm1/(4*dx**2) - qajm2/(16*dx**2) + 3*qajp1/(4*dx**2) - qajp2/(16*dx**2)
    d = -17*qajm1/(24*dx) + 5*qajm2/(48*dx) + 17*qajp1/(24*dx) - 5*qajp2/(48*dx)
    e = 1067*qaj/960 - 29*qajm1/480 + 3*qajm2/640 - 29*qajp1/480 + 3*qajp2/640

    
    SI = 77051*qaj**2/1680 - 24923*qaj*qajm1/420 + 7547*qaj*qajm2/560 - 24923*qaj*qajp1/420 + 7547*qaj*qajp2/560 + 104963*qajm1**2/5040 - 51001*qajm1*qajm2/5040 + 89549*qajm1*qajp1/2520 - 38947*qajm1*qajp2/5040 + 1727*qajm2**2/1260 - 38947*qajm2*qajp1/5040 + 8209*qajm2*qajp2/5040 + 104963*qajp1**2/5040 - 51001*qajp1*qajp2/5040 + 1727*qajp2**2/1260
    
    return a,b,c,d,e,SI

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


##### Reconstruction Functions #####
"""
This group of functions takes the values of the cell average quantity in the relevant cells and produces the coefficients of the cubic that approximates the function in the jth cell

The Recon function then uses the basic coefficeint recontruction function ReconCoeffs to reconstruct x_{j+1/2} and x_{j-1/2} for each cell in the domain.
With constant values used in the ghost cells
"""
def DegreeP(q,eps):
    n = len(q)
    m = n
    for i in range(n):
        if( abs(q[i]) < eps):
            m = m -1
        else:
            break
    return m

def SAL(q):
    qmin = min(q)
    qmax = max(q)
    if qmin>0:
        return 1
    if qmax<0 :
        return 1
    else:
        return 0

def minmodL(q):
    qmin = min(q)
    qmax = max(q)
    if qmin>0:
        return qmin
    if qmax<0 :
        return qmax
    else:
        return 0


def ReconCoeffs(qajm4,qajm3,qajm2,qajm1,qaj,qajp1,qajp2,qajp3,qajp4,dx):

    eps = 10.0**(-12)
    #Start with constant
    #Cubic
    pjm3toja33,pjm3toja32,pjm3toja31,pjm3toja30,Bjm3toj = P3jm3toj(qaj,qajm1,qajm2,qajm3,dx)
    pjm2tojp1a33,pjm2tojp1a32,pjm2tojp1a31,pjm2tojp1a30,Bjm2tojp1 = P3jm2tojp1(qaj,qajm1,qajm2,qajp1,dx)
    pjm1tojp2a33,pjm1tojp2a32,pjm1tojp2a31,pjm1tojp2a30,Bjm1tojp2 = P3jm1tojp2(qaj,qajm1,qajp1,qajp2,dx)
    pjtojp3a33,pjtojp3a32,pjtojp3a31,pjtojp3a30,Bjtojp3 = P3jtojp3(qaj,qajp1,qajp2,qajp3,dx)

    
    #ExpLeftIn
    pjm4tojm1a33,pjm4tojm1a32,pjm4tojm1a31,pjm4tojm1a30,pjm4tojm1SI = P3jm3toj(qajm1,qajm2,qajm3,qajm4,dx)

    #ExpRightIn
    pjp1tojp4a33,pjp1tojp4a32,pjp1tojp4a31,pjp1tojp4a30,pjp1tojp4SI = P3jtojp3(qajp1,qajp2,qajp3,qajp4,dx)
    
    
    RIxjph = pjp1tojp4a33/4*(-dx/2)**4 + pjp1tojp4a32/3*(-dx/2)**3 + pjp1tojp4a31/2*(-dx/2)**2 + pjp1tojp4a30*(-dx/2)
    RIxjmh = pjp1tojp4a33/4*(-3*dx/2)**4 + pjp1tojp4a32/3*(-3*dx/2)**3 + pjp1tojp4a31/2*(-3*dx/2)**2 + pjp1tojp4a30*(-3*dx/2) 
    PlusCAErr = abs((RIxjph -  RIxjmh)/dx-qaj )


    RIxjph = pjm4tojm1a33/4*(3*dx/2)**4 + pjm4tojm1a32/3*(3*dx/2)**3 + pjm4tojm1a31/2*(3*dx/2)**2 + pjm4tojm1a30*(3*dx/2)
    RIxjmh = pjm4tojm1a33/4*(dx/2)**4 + pjm4tojm1a32/3*(dx/2)**3 + pjm4tojm1a31/2*(dx/2)**2 + pjm4tojm1a30*(dx/2) 
    MinusCAErr = abs((RIxjph -  RIxjmh)/dx-qaj )
    
    
    P3jm3tojD = DegreeP([pjm3toja33,pjm3toja32,pjm3toja31,pjm3toja30],eps)
    P3jm2tojp1D = DegreeP([pjm2tojp1a33,pjm2tojp1a32,pjm2tojp1a31,pjm2tojp1a30],eps)
    P3jm1tojp2D =DegreeP([pjm1tojp2a33,pjm1tojp2a32,pjm1tojp2a31,pjm1tojp2a30],eps)
    P3jtojp3D =  DegreeP([pjtojp3a33,pjtojp3a32,pjtojp3a31,pjtojp3a30],eps)
    
    if (abs(PlusCAErr -MinusCAErr)<eps):
        
        #Pick lowest order one
        if( min(P3jm3tojD,P3jm2tojp1D,P3jm1tojp2D,P3jtojp3D)< 4):
            if(P3jm3tojD <= min(P3jm2tojp1D,P3jm1tojp2D,P3jtojp3D) ):
                Ra3 = pjm3toja33
                Ra2 = pjm3toja32
                Ra1 = pjm3toja31
                Ra0 = pjm3toja30
            elif(P3jm2tojp1D <= min(P3jm3tojD,P3jm1tojp2D,P3jtojp3D) ):
                Ra3 = pjm2tojp1a33
                Ra2 = pjm2tojp1a32
                Ra1 = pjm2tojp1a31
                Ra0 = pjm2tojp1a30
            elif(P3jm1tojp2D <= min(P3jm3tojD,P3jm2tojp1D,P3jtojp3D) ):
                Ra3 = pjm1tojp2a33
                Ra2 = pjm1tojp2a32
                Ra1 = pjm1tojp2a31
                Ra0 = pjm1tojp2a30
            else:
                Ra3 = pjm3toja33
                Ra2 = pjm3toja32
                Ra1 = pjm3toja31
                Ra0 = pjm3toja30
        else:     
            a66,a65,a64,Ra3,Ra2,Ra1,Ra0,SI = P6jm3tojp3(qajm3,qajm2,qajm1,qaj,qajp1,qajp2,qajp3,dx)
    elif(PlusCAErr > MinusCAErr):
        Ra3 = pjm3toja33
        Ra2 = pjm3toja32
        Ra1 = pjm3toja31
        Ra0 = pjm3toja30
    else:
        Ra3 = pjtojp3a33
        Ra2 = pjtojp3a32
        Ra1 = pjtojp3a31
        Ra0 = pjtojp3a30
    
        
    return 0,Ra3,Ra2,Ra1,Ra0

def Recon(qA,x,dx,nGcells,eps):
    n = len(x)
    qjphs = zeros(n)
    xjphs = zeros(n)
    qjmhs = zeros(n)
    xjmhs = zeros(n)
    
    for i in range(nGcells,n-nGcells):
        a4,a3,a2,a1,a0= ReconCoeffs(qA[i-4],qA[i-3],qA[i-2],qA[i-1],qA[i],qA[i+1],qA[i+2],qA[i+3],qA[i+4],dx)

        
        #Reconstruction at  x_{j+1/2}
        qjphs[i] = a4*(dx/2)**4 + a3*(dx/2)**3 + a2*(dx/2)**2 + a1*(dx/2) + a0
        xjphs[i] = x[i] + 0.5*dx
        
        #Reconstruction at  x_{j-1/2}
        qjmhs[i] = a4*(-dx/2)**4 +a3*(-dx/2)**3 + a2*(-dx/2)**2 + a1*(-dx/2) + a0
        xjmhs[i] = x[i] - 0.5*dx
    
    for i in range(nGcells):
        qjphs[i] = qA[i]
        xjphs[i] = x[i] + 0.5*dx
        qjmhs[i] = qA[i]
        xjmhs[i] = x[i] - 0.5*dx

        qjphs[n- 1 -i] = qA[n- 1 -i]
        xjphs[n- 1 -i] = x[n- 1 -i] + 0.5*dx
        qjmhs[n- 1 -i] = qA[n- 1 -i]
        xjmhs[n- 1 -i] = x[n- 1 -i] - 0.5*dx
        
    return xjphs,qjphs,xjmhs,qjmhs


##### Initial Condition Functions #####
"""
PolynomialEval - given a list of coefficients evaluates the polynomial with thos coefficients at x

IntPolynomialEval- given a list of coefficients evaluates the antiderivative of the polynomial with thos coefficients at x

PP_IC - given the cell centers x, cell width dx, and the Coefficients for the left polynomial (for when x < 0) and the coefficients of the right polynomial (for when x => 0)
"""

def PolynomialEval(Coeffs,x):
    n = len(Coeffs)
    Sum = Coeffs[0]
    for i in range(1,n):
        Sum = Sum + Coeffs[i]*x**i
    return Sum

def IntPolynomialEval(Coeffs,x):
    n = len(Coeffs)
    Sum = Coeffs[0]*x
    for i in range(1,n):
        Sum = Sum + Coeffs[i]/(i+1)*x**(i+1)
    return Sum


def PP_IC(x,dx,LeftCoeffs,RightCoeffs):
    n= len(x)
    hA = zeros(n)
    hMjph = zeros(n)
    hMjmh = zeros(n)
    xjphs = zeros(n)
    xjmhs = zeros(n)
    
    for i in range(n):
        xiph = x[i] + 0.5*dx
        xi = x[i]
        ximh = x[i] - 0.5*dx
        xjphs[i] = xiph
        xjmhs[i] = ximh
        if (xiph > 0 and ximh < 0):
            hMjph[i] = PolynomialEval(RightCoeffs,xiph)
            hMjmh[i] = PolynomialEval(LeftCoeffs,ximh)
            
            hIph = IntPolynomialEval(LeftCoeffs,xi)
            hImh = IntPolynomialEval(LeftCoeffs,ximh)
            hALeft = (hIph - hImh)/(0.5*dx)
            
            hIph = IntPolynomialEval(RightCoeffs,xiph)
            hImh = IntPolynomialEval(RightCoeffs,xi)
            hARight = (hIph - hImh)/(0.5*dx)
            hA[i] = hALeft + hARight

            
        elif(xiph <0):
            
            hMjph[i] = PolynomialEval(LeftCoeffs,xiph)
            hMjmh[i] = PolynomialEval(LeftCoeffs,ximh)
            
            hIph = IntPolynomialEval(LeftCoeffs,xiph)
            hImh = IntPolynomialEval(LeftCoeffs,ximh)
            hA[i] = (hIph - hImh)/dx

        else:
            hMjph[i] = PolynomialEval(RightCoeffs,xiph)
            hMjmh[i] = PolynomialEval(RightCoeffs,ximh)
            
            hIph = IntPolynomialEval(RightCoeffs,xiph)
            hImh = IntPolynomialEval(RightCoeffs,ximh)
            hA[i] = (hIph - hImh)/dx
            
    return hA,hMjmh,hMjph,xjmhs,xjphs


#Main program
g = 9.81

a0 = 10
a1 = 2.0/8.0
a2 = 1.0/20.0
a3 = -19.0/400
a4 = -3.0/1000
a5 =147.0/80000
a6 =9.0/200000
a7 =-3.0/125000
a8 =-1.0/5000000
a9 =1.0/10000000

# a0 = 10
# a1 = 2.0/8.0
# a2 = 0
# a3 = 0
# a4 = 0
# a5 =0
# a6 =0
# a7 =0
# a8 =0
# a9 =0

#Coefficients for left and right polynomials (can make discontinuity, by negating coefficients in right polynomial)
LeftCoeffs = [a0,a1,a2,a3,a4,a5,a6,a7,a8,a9]
RightCoeffs = [-a0,-a1,-a2,-a3,-a4,-a5,-a6,-a7,-a8,-a9]

# LeftCoeffs = [a0,a1,a2,a3,a4,a5]
# RightCoeffs = [a0,a1,a2,a3,a4,a5]

eps = 10.0**(-10)
nGcells = 5
expn = 10
lown = 40

dxs = zeros(expn)
L2phs = zeros(expn)
L2mhs  = zeros(expn)

#Loops for a variety of dx values
for expi in range(expn):
    ncurr = lown*(2**expi)
    
    sx = -10.0
    ex = 10.0
    dx = ((ex - sx))/(ncurr-1)
    x = arange(sx ,ex+(0.1)*dx,dx)
    nx = len(x)
    
    hA,hMjmha,hMjpha,xjmhsa,xjphsa = PP_IC(x,dx,LeftCoeffs,RightCoeffs)
    xjphs,hjphs,xjmhs,hjmhs = Recon(hA,x,dx,nGcells,eps)
    
    dxs[expi] = dx
    L2phs[expi] = norm(hjphs[nGcells:-nGcells] - hMjpha[nGcells:-nGcells],ord=2 )/ norm(hMjpha[nGcells:-nGcells],ord=2 )
    L2mhs[expi] = norm(hjmhs[nGcells:-nGcells] - hMjmha[nGcells:-nGcells],ord=2 )/ norm(hMjmha[nGcells:-nGcells],ord=2 )

#Plot the convergence plots
loglog(dxs,L2phs,'ob',label='Recon at x^{-}_{j+1/2}')
loglog(dxs,L2mhs,'+r',label='Recon at x^{+}_{j-1/2}')
loglog(dxs,(L2phs[0]/dxs[0]**1)*array(dxs)**1,'-r',label='Slope - 1')
loglog(dxs,(L2phs[0]/dxs[0]**2)*array(dxs)**2,'-b',label='Slope - 2')
loglog(dxs,(L2phs[0]/dxs[0]**3)*array(dxs)**3,'-g',label='Slope - 3')
loglog(dxs,(L2phs[0]/dxs[0]**4)*array(dxs)**4,'-k',label='Slope - 4')

legend()
s = "Convergence \n\n"
s = s + "Left Coefficients : " + str(LeftCoeffs ) + "\n"
s = s + "Right Coefficients : " + str(RightCoeffs ) + "\n"
title(s)
ylabel('L2')
xlabel('dx')