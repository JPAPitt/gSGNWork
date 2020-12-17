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


##### Reconstruction Functions #####
"""
This group of functions takes the values of the cell average quantity in the relevant cells and produces the coefficients of the cubic that approximates the function in the jth cell

The Recon function then uses the basic coefficeint recontruction function ReconCoeffs to reconstruct x_{j+1/2} and x_{j-1/2} for each cell in the domain.
With constant values used in the ghost cells
"""

def minmodL(q):
    qmin = min(q)
    qmax = max(q)
    if qmin>0:
        return qmin
    if qmax<0 :
        return qmax
    else:
        return 0


def ReconCoeffs(qajm3,qajm2,qajm1,qaj,qajp1,qajp2,qajp3,dx):
    #Start with constant
    P0ja0 = qaj

    #Linear
    P1jtojp1a11,P1jtojp1a10,P1jtojp1SI  =  P1jtojp1(qaj,qajp1,dx)
    P1jm1toja11,P1jm1toja10,P1jm1tojSI  =  P1jm1toj(qaj,qajm1,dx)
            
  
    #Quadratic
    P2jtojp2a22,P2jtojp2a21,P2jtojp2a20,P2jtojp2SI =   P2jtojp2(qaj,qajp1,qajp2,dx)
    P2jm1tojp1a22,P2jm1tojp1a21,P2jm1tojp1a20,P2jm1tojp1SI =   P2jm1tojp1(qajp1,qaj,qajm1,dx)
    P2jm2toja22,P2jm2toja21,P2jm2toja20,P2jm2tojSI = P2jm2toj(qajm2,qajm1,qaj,dx)   
    
    #Cubic
    pjm3toja,pjm3tojb,pjm3tojc,pjm3tojd,Bjm3toj = P3jm3toj(qaj,qajm1,qajm2,qajm3,dx)
    pjm2tojp1a,pjm2tojp1b,pjm2tojp1c,pjm2tojp1d,Bjm2tojp1 = P3jm2tojp1(qaj,qajm1,qajm2,qajp1,dx)
    pjm1tojp2a,pjm1tojp2b,pjm1tojp2c,pjm1tojp2d,Bjm1tojp2 = P3jm1tojp2(qaj,qajm1,qajp1,qajp2,dx)
    pjtojp3a,pjtojp3b,pjtojp3c,pjtojp3d,Bjtojp3 = P3jtojp3(qaj,qajp1,qajp2,qajp3,dx)
    
    TotalSI =  P1jtojp1SI + P1jm1tojSI
    if TotalSI > 10.0**(-10):
        P1Ra1 = (P1jtojp1SI/TotalSI)*P1jm1toja11 + (P1jm1tojSI/TotalSI)*P1jtojp1a11 
        P1Ra0 = (P1jtojp1SI/TotalSI)*P1jm1toja10 + (P1jm1tojSI/TotalSI)*P1jtojp1a10 
    else:
        P1Ra1 = 0
        P1Ra0 = qaj

     
        
    TotalSI =  (P2jtojp2SI + P2jm1tojp1SI + P2jm2tojSI)*2.0
    if TotalSI > 10.0**(-10):
        P2Ra2 = (((P2jm1tojp1SI + P2jtojp2SI) /TotalSI)*P2jm2toja22 + ((P2jm2tojSI + P2jtojp2SI) /TotalSI)*P2jm1tojp1a22 +  ((P2jm2tojSI + P2jm1tojp1SI) /TotalSI)*P2jtojp2a22)
        P2Ra1 = (((P2jm1tojp1SI + P2jtojp2SI)/TotalSI)*P2jm2toja21 + ((P2jm2tojSI + P2jtojp2SI) /TotalSI)*P2jm1tojp1a21 +  ((P2jm2tojSI + P2jm1tojp1SI) /TotalSI)*P2jtojp2a21)
        P2Ra0 = (((P2jm1tojp1SI + P2jtojp2SI) /TotalSI)*P2jm2toja20 + ((P2jm2tojSI + P2jtojp2SI) /TotalSI)*P2jm1tojp1a20 +  ((P2jm2tojSI + P2jm1tojp1SI)/TotalSI)*P2jtojp2a20 )
    else:
        P2Ra2 = 0
        P2Ra1 = 0
        P2Ra0 = qaj           
    
    
    TotalSI =  (Bjm3toj + Bjm2tojp1  +  Bjm1tojp2  +Bjtojp3 )*3.0
    if TotalSI > 10.0**(-14):
        P3Ra3 = ((Bjm2tojp1  +  Bjm1tojp2  +Bjtojp3) /TotalSI)*pjm3toja \
                + ((Bjm3toj  +  Bjm1tojp2  +Bjtojp3) /TotalSI)*pjm2tojp1a \
                + ((Bjm3toj  +  Bjm2tojp1  +Bjtojp3) /TotalSI)*pjm1tojp2a \
                + ((Bjm3toj  +  Bjm2tojp1  +Bjm1tojp2) /TotalSI)*pjtojp3a
                
        P3Ra2 = ((Bjm2tojp1  +  Bjm1tojp2  +Bjtojp3) /TotalSI)*pjm3tojb \
                + ((Bjm3toj  +  Bjm1tojp2  +Bjtojp3) /TotalSI)*pjm2tojp1b \
                + ((Bjm3toj  +  Bjm2tojp1  +Bjtojp3) /TotalSI)*pjm1tojp2b \
                + ((Bjm3toj  +  Bjm2tojp1  +Bjm1tojp2) /TotalSI)*pjtojp3b                     
                 
        P3Ra1 = ((Bjm2tojp1  +  Bjm1tojp2  +Bjtojp3) /TotalSI)*pjm3tojc \
                + ((Bjm3toj  +  Bjm1tojp2  +Bjtojp3) /TotalSI)*pjm2tojp1c \
                + ((Bjm3toj  +  Bjm2tojp1  +Bjtojp3) /TotalSI)*pjm1tojp2c \
                + ((Bjm3toj  +  Bjm2tojp1  +Bjm1tojp2) /TotalSI)*pjtojp3c

        P3Ra0 = ((Bjm2tojp1  +  Bjm1tojp2  +Bjtojp3) /TotalSI)*pjm3tojd \
                + ((Bjm3toj  +  Bjm1tojp2  +Bjtojp3) /TotalSI)*pjm2tojp1d \
                + ((Bjm3toj  +  Bjm2tojp1  +Bjtojp3) /TotalSI)*pjm1tojp2d \
                + ((Bjm3toj  +  Bjm2tojp1  +Bjm1tojp2) /TotalSI)*pjtojp3d
                
    else:
        P3Ra3 = 0
        P3Ra2 = 0
        P3Ra1 = 0
        P3Ra0 = qaj  
    
    
    if (P1Ra1 == 0):
        nRa3 = 0
        nRa2 = 0
        nRa1 = 0
        nRa0 = qaj  
    elif(P2Ra2 == 0):
        nRa3 = 0
        nRa2 = 0
        nRa1 = P1Ra1
        nRa0 = P1Ra0 
    elif(P3Ra3 == 0):
        nRa3 = 0
        nRa2 = P2Ra2
        nRa1 = P2Ra1
        nRa0 = P2Ra0 
    else:
        nRa3 = P3Ra3
        nRa2 = P3Ra2
        nRa1 = P3Ra1
        nRa0 = P3Ra0 
    
    return P3Ra3,P3Ra2,P3Ra1,P3Ra0

def Recon(qA,x,dx,nGcells,eps):
    n = len(x)
    qjphs = zeros(n)
    xjphs = zeros(n)
    qjmhs = zeros(n)
    xjmhs = zeros(n)
    
    for i in range(nGcells,n-nGcells):
        a3,a2,a1,a0 = ReconCoeffs(qA[i-3],qA[i-2],qA[i-1],qA[i],qA[i+1],qA[i+2],qA[i+3],dx)

        
        #Reconstruction at  x_{j+1/2}
        qjphs[i] = a3*(dx/2)**3 + a2*(dx/2)**2 + a1*(dx/2) + a0
        xjphs[i] = x[i] + 0.5*dx
        
        #Reconstruction at  x_{j-1/2}
        qjmhs[i] = a3*(-dx/2)**3 + a2*(-dx/2)**2 + a1*(-dx/2) + a0
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
        if(xiph <0):
            
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

#Coefficients for left and right polynomials (can make discontinuity, by negating coefficients in right polynomial)
LeftCoeffs = [a0,a1,a2,a3,a4,a5,a6,a7,a8,a9]
RightCoeffs = [a0,a1,a2,a3,a4,a5,a6,a7,a8,a9]

eps = 10.0**(-10)
nGcells = 4
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
    x = arange(sx+0.1*dx ,ex+(0.1)*dx,dx)
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