"""
This is a python script that tests reconstructions - in particular this one is the WENO with cubic base stencils which at most 7th order accurate in smooth regions (due to using)
weights to match the polynomial that covers all cells used in the candidiate stencils from j-3,....,j+3. These weights are unique to each point in the cell, hence the need to
have different optimal weights for j+1/2 and j-1/2 in a cell, also if we were to approximate derivaitves optimal weights for the derivative approximations at points would
also be required.

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
This group of functions takes the values of the cell average quantity in the relevant cells and produces the reconstruction at x_{j+1/2} for Reconjph and x_{j-1/2} for Reconjmh 

The Recon function then uses the two basic WENO recontruction functions Reconjph and Reconjmh to reconstruct x_{j+1/2} and x_{j-1/2} for each cell in the domain.
With constant values used in the ghost cells
"""
def Reconjph(qajm3,qajm2,qajm1,qaj,qajp1,qajp2,qajp3,dx,eps):

    pjm3toja,pjm3tojb,pjm3tojc,pjm3tojd,Bjm3toj = P3jm3toj(qaj,qajm1,qajm2,qajm3,dx)
    pjm2tojp1a,pjm2tojp1b,pjm2tojp1c,pjm2tojp1d,Bjm2tojp1 = P3jm2tojp1(qaj,qajm1,qajm2,qajp1,dx)
    pjm1tojp2a,pjm1tojp2b,pjm1tojp2c,pjm1tojp2d,Bjm1tojp2 = P3jm1tojp2(qaj,qajm1,qajp1,qajp2,dx)
    pjtojp3a,pjtojp3b,pjtojp3c,pjtojp3d,Bjtojp3 = P3jtojp3(qaj,qajp1,qajp2,qajp3,dx)
    
            
    iw1 = ((1.0/35.0) / (eps +Bjm3toj )**2)
    iw2 = ((12.0/35.0) / (eps +Bjm2tojp1 )**2)
    iw3 = ((18.0/35.0)  / (eps +Bjm1tojp2 )**2)
    iw4 = ((4.0/35.0)  / (eps +Bjtojp3 )**2)

    
    w1 = iw1 / (iw1 + iw2 + iw3 + iw4 )
    w2 = iw2 / (iw1 + iw2 + iw3 + iw4 )
    w3 = iw3 / (iw1 + iw2 + iw3 + iw4 )
    w4 = iw4 / (iw1 + iw2 + iw3 + iw4 )
    
    qa = w1*pjm3toja + w2*pjm2tojp1a + w3*pjm1tojp2a + w4*pjtojp3a
    qb = w1*pjm3tojb + w2*pjm2tojp1b + w3*pjm1tojp2b + w4*pjtojp3b
    qc = w1*pjm3tojc + w2*pjm2tojp1c + w3*pjm1tojp2c + w4*pjtojp3c
    qd = w1*pjm3tojd + w2*pjm2tojp1d + w3*pjm1tojp2d + w4*pjtojp3d
        

    return qa*(dx/2)**3 + qb*(dx/2)**2 + qc*(dx/2) + qd


def Reconjmh(qajm3,qajm2,qajm1,qaj,qajp1,qajp2,qajp3,dx,eps):

    pjm3toja,pjm3tojb,pjm3tojc,pjm3tojd,Bjm3toj = P3jm3toj(qaj,qajm1,qajm2,qajm3,dx)
    pjm2tojp1a,pjm2tojp1b,pjm2tojp1c,pjm2tojp1d,Bjm2tojp1 = P3jm2tojp1(qaj,qajm1,qajm2,qajp1,dx)
    pjm1tojp2a,pjm1tojp2b,pjm1tojp2c,pjm1tojp2d,Bjm1tojp2 = P3jm1tojp2(qaj,qajm1,qajp1,qajp2,dx)
    pjtojp3a,pjtojp3b,pjtojp3c,pjtojp3d,Bjtojp3 = P3jtojp3(qaj,qajp1,qajp2,qajp3,dx)
    
    iw1 = ((4.0/35.0) / (eps +Bjm3toj )**2)
    iw2 = ((18.0/35.0) / (eps +Bjm2tojp1 )**2)
    iw3 = ((12.0/35.0)  / (eps +Bjm1tojp2 )**2)
    iw4 = ((1.0/35.0)  / (eps +Bjtojp3 )**2)
    
    
    w1 = iw1 / (iw1 + iw2 + iw3 + iw4 )
    w2 = iw2 / (iw1 + iw2 + iw3 + iw4 )
    w3 = iw3 / (iw1 + iw2 + iw3 + iw4 )
    w4 = iw4 / (iw1 + iw2 + iw3 + iw4 )
    
    qa = w1*pjm3toja + w2*pjm2tojp1a + w3*pjm1tojp2a + w4*pjtojp3a
    qb = w1*pjm3tojb + w2*pjm2tojp1b + w3*pjm1tojp2b + w4*pjtojp3b
    qc = w1*pjm3tojc + w2*pjm2tojp1c + w3*pjm1tojp2c + w4*pjtojp3c
    qd = w1*pjm3tojd + w2*pjm2tojp1d + w3*pjm1tojp2d + w4*pjtojp3d
        

    return qa*(-dx/2)**3 + qb*(-dx/2)**2 + qc*(-dx/2) + qd

def Recon(qA,x,dx,nGcells,eps):
    n = len(x)
    qjphs = zeros(n)
    xjphs = zeros(n)
    qjmhs = zeros(n)
    xjmhs = zeros(n)
    
    for i in range(nGcells,n-nGcells):
        
        #Reconstruction at  x_{j+1/2}
        qjphs[i] = Reconjph(qA[i-3],qA[i-2],qA[i-1],qA[i],qA[i+1],qA[i+2],qA[i+3],dx,eps)
        xjphs[i] = x[i] + 0.5*dx
        
        #Reconstruction at  x_{j-1/2}
        qjmhs[i] = Reconjmh(qA[i-3],qA[i-2],qA[i-1],qA[i],qA[i+1],qA[i+2],qA[i+3],dx,eps)
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
RightCoeffs = [-a0,-a1,-a2,-a3,a4,a5,a6,a7,a8,a9]

eps = 10.0**(-12)
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
loglog(dxs,(min(L2phs[0],L2mhs[0])/dxs[0]**3)*array(dxs)**3,'-r',label='Slope - 3')
loglog(dxs,(min(L2phs[0],L2mhs[0])/dxs[0]**4)*array(dxs)**4,'-b',label='Slope - 4')
loglog(dxs,(min(L2phs[0],L2mhs[0])/dxs[0]**5)*array(dxs)**5,'-g',label='Slope - 5')
loglog(dxs,(min(L2phs[0],L2mhs[0])/dxs[0]**6)*array(dxs)**6,'-k',label='Slope - 6')
loglog(dxs,(min(L2phs[0],L2mhs[0])/dxs[0]**7)*array(dxs)**7,'-r',label='Slope - 7')

legend()
s = "Convergence \n\n"
s = s + "Left Coefficients : " + str(LeftCoeffs ) + "\n"
s = s + "Right Coefficients : " + str(RightCoeffs ) + "\n"
title(s)
ylabel('L2')
xlabel('dx')