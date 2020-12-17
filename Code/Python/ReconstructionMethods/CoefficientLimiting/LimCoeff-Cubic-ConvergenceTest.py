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
    return (qajp1 - qaj)/ dx ,qaj

def P1jm1toj(qaj,qajm1,dx):
    return (qaj - qajm1)/ dx ,qaj

def P1jWeird(qajp1,qaj,qajm1,dx):
    return (qajp1 - qajm1)/ (2*dx) ,qaj

def P2jtojp2(qaj,qajp1,qajp2,dx):
    a = (qajp2 - 2*qajp1 + qaj)/ (2*dx**2)
    b = (-qajp2 + 4*qajp1 - 3*qaj)/ (2*dx)
    c =  -qajp2/24 + qajp1/12 + 23*qaj/24
    return a,b,c

def P2jm1tojp1(qajm1,qaj,qajp1,dx):
    a = (qajp1 - 2*qaj + qajm1)/ (2*dx**2)
    b = (qajp1 - qajm1)/ (2*dx)
    c =  13*qaj/12 - qajp1/24 - qajm1/24
    return a,b,c

def P2jm2toj(qajm2,qajm1,qaj,dx):
    a = (qaj - 2*qajm1 + qajm2)/ (2*dx**2)
    b = (qajm2 - 4*qajm1 + 3*qaj)/ (2*dx)
    c =  -qajm2/24 + qajm1/12 + 23*qaj/24
    return a,b,c


def P3jtojp3(qaj,qajp1,qajp2,qajp3,dx):
    a=(3*qajp1 - 3*qajp2 + qajp3 - qaj)/(6*dx**3)
    b=(-5*qajp1 + 4*qajp2 - qajp3 + 2*qaj)/(2*dx**2)
    c=(69*qajp1 - 33*qajp2 + 7*qajp3 - 43*qaj)/(24*dx)
    d=5*qajp1/24 - qajp2/6 + qajp3/24 + 11*qaj/12
    return a,b,c,d


def P3jm1tojp2(qaj,qajm1,qajp1,qajp2,dx):
    a=(-3*qajp1 + qajp2 - qajm1 + 3*qaj)/(6*dx**3)
    b=(qajp1 + qajm1 - 2*qaj)/(2*dx**2)
    c=(27*qajp1 - 5*qajp2 - 7*qajm1 - 15*qaj)/(24*dx)
    d=-qajp1/24 - qajm1/24 + 13*qaj/12
    return a,b,c,d

def P3jm2tojp1(qaj,qajm1,qajm2,qajp1,dx):
    a=(qajp1 + 3*qajm1 - qajm2 - 3*qaj)/(6*dx**3)
    b=(qajp1 + qajm1 - 2*qaj)/(2*dx**2)
    c=(7*qajp1- 27*qajm1 + 5*qajm2 + 15*qaj)/(24*dx)
    d=-qajp1/24 - qajm1/24 + 13*qaj/12
    return a,b,c,d


def P3jm3toj(qaj,qajm1,qajm2,qajm3,dx):
    a=(-3*qajm1 + 3*qajm2 - qajm3 + qaj)/(6*dx**3)
    b=(-5*qajm1 + 4*qajm2 - qajm3 + 2*qaj)/(2*dx**2)
    c=(-69*qajm1 + 33*qajm2 - 7*qajm3 + 43*qaj)/(24*dx)
    d=5*qajm1/24 - qajm2/6 + qajm3/24 + 11*qaj/12
    return a,b,c,d


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
    P3r0a33,P3r0a32,P3r0a31,P3r0a30 = P3jtojp3(qaj,qajp1,qajp2,qajp3,dx)
    P3r1a33,P3r1a32,P3r1a31,P3r1a30 = P3jm1tojp2(qaj,qajm1,qajp1,qajp2,dx)
    P3r2a33,P3r2a32,P3r2a31,P3r2a30 = P3jm2tojp1(qaj,qajm1,qajm2,qajp1,dx)
    P3r3a33,P3r3a32,P3r3a31,P3r3a30 =  P3jm3toj(qaj,qajm1,qajm2,qajm3,dx)
    
    P2r0a22,P2r0a21,P2r0a20=   P2jtojp2(qaj,qajp1,qajp2,dx)
    P2r1a22,P2r1a21,P2r1a20  =  P2jm1tojp1(qajm1,qaj,qajp1,dx)
    P2r2a22,P2r2a21,P2r2a20  =  P2jm2toj(qajm2,qajm1,qaj,dx)
    
    P1r0a11,P1r0a10  =  P1jtojp1(qaj,qajp1,dx)
    P1r1a11,P1r1a10  =  P1jm1toj(qaj,qajm1,dx)
    
    P1Wa11,P1Wa10  = P1jWeird(qajp1,qaj,qajm1,dx)
    
    # w1 = minmodL([P1r0a11,P1r1a11,P1Wa11])  / P1Wa11
    # P1a11 = P1Wa11
    P1a11w1 = minmodL([P1r0a11,P1r1a11,P1Wa11])
    
    #w2 = minmodL([P2r0a22,P2r1a2,P2r2a22])  / P2r1a22
    #P2a22 = P2r1a22
    P2a22w2 = minmodL([P2r0a22,P2r1a22,P2r2a22]) 
    
    #w3 = minmodL([P3r0a33,P3r1a33,P3r2a33,P3r3a33])  / P3r2a33
    #P3a33 = P3r2a33
    P3a33w3 = minmodL([P3r0a33,P3r1a33,P3r2a33,P3r3a33])
    
    #Weight on second coefficient
    if ( abs(P3r2a33) > 10.0**(-12.0)):
        w3 = minmodL([P3r0a33,P3r1a33,P3r2a33,P3r3a33])  / P3r2a33
    else:
        w3 = 0
    #P3a32w3 = w3*P3r2a32
    P3a32w3Corr = w3*(P3r2a32 - P2a22w2) + P2a22w2
    
    
    #Weight on first coefficient, want it to reduce to Reconed P2a22
    if (abs(P2r1a22)  > 10.0**(-12.0)):
        w2 = minmodL([P2r0a22,P2r1a22,P2r2a22])  / P2r1a22
    else:
        w2 = 0
    P2a21w2 = w2 *P2r1a21
    
    P2a21w2Corr = w2*(P2r1a21 - P1a11w1) + P1a11w1
    
    P3a31w3Corr = w3*(P3r2a31 - P2a21w2Corr) + P2a21w2Corr
    
    na33 = P3a33w3
    na32 = P3a32w3Corr
    na31 = P3a31w3Corr
    na30 = qaj - na32/3*(dx/2)**2 
    
    return na33,na32,na31,na30

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