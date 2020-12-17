"""
This is a python script that tests reconstructions - in particular this one chooses the reconstruction using a minmod comparison of the coefficients of the possible
forward and backward polynomials as well as the higher order centered one using all the cells.

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

def P2jm1tojp1(qajp1,qaj,qajm1,dx):
    a =-qaj/dx**2 + qajm1/(2*dx**2) + qajp1/(2*dx**2)
    b=-qajm1/(2*dx) + qajp1/(2*dx)
    c=13*qaj/12 - qajm1/24 - qajp1/24
    return a,b,c

#
def P2jtojp2(qaj,qajp1,qajp2,dx):
    a = (qajp2 - 2*qajp1 + qaj)/ (2*dx**2)
    b = (-qajp2 + 4*qajp1 - 3*qaj)/ (2*dx)
    c =  -qajp2/24 + qajp1/12 + 23*qaj/24
    return a,b,c

def P2jm2toj(qajm2,qajm1,qaj,dx):
    a = (qaj - 2*qajm1 + qajm2)/ (2*dx**2)
    b = (qajm2 - 4*qajm1 + 3*qaj)/ (2*dx)
    c =  -qajm2/24 + qajm1/12 + 23*qaj/24
    return a,b,c

def P4jm2tojp2(qaj,qajm1,qajm2,qajp1,qajp2,dx):
    a =-qaj/(4*dx**4) - qajm1/(6*dx**4) + qajm2/(24*dx**4) - qajp1/(6*dx**4) + qajp2/(24*dx**4)
    b= qajm1/(6*dx**3) - qajm2/(12*dx**3) - qajp1/(6*dx**3) + qajp2/(12*dx**3)
    c = -11*qaj/(8*dx**2) + 3*qajm1/(4*dx**2) - qajm2/(16*dx**2) + 3*qajp1/(4*dx**2) - qajp2/(16*dx**2)
    d = -17*qajm1/(24*dx) + 5*qajm2/(48*dx) + 17*qajp1/(24*dx) - 5*qajp2/(48*dx)
    e = 1067*qaj/960 - 29*qajm1/480 + 3*qajm2/640 - 29*qajp1/480 + 3*qajp2/640
    return a,b,c,d,e



def P3jtojp3(qaj,qajp1,qajp2,qajp3,dx):
    a=(3*qajp1 - 3*qajp2 + qajp3 - qaj)/(6*dx**3)
    b=(-5*qajp1 + 4*qajp2 - qajp3 + 2*qaj)/(2*dx**2)
    c=(69*qajp1 - 33*qajp2 + 7*qajp3 - 43*qaj)/(24*dx)
    d=5*qajp1/24 - qajp2/6 + qajp3/24 + 11*qaj/12
    return a,b,c,d

def P3jm3toj(qaj,qajm1,qajm2,qajm3,dx):
    a=(-3*qajm1 + 3*qajm2 - qajm3 + qaj)/(6*dx**3)
    b=(-5*qajm1 + 4*qajm2 - qajm3 + 2*qaj)/(2*dx**2)
    c=(-69*qajm1 + 33*qajm2 - 7*qajm3 + 43*qaj)/(24*dx)
    d=5*qajm1/24 - qajm2/6 + qajm3/24 + 11*qaj/12
    return a,b,c,d

def P6jm3tojp3(qaj,qajm1,qajm2,qajm3,qajp1,qajp2,qajp3,dx):
    a =-qaj/(36*dx**6) + qajm1/(48*dx**6) - qajm2/(120*dx**6) + qajm3/(720*dx**6) + qajp1/(48*dx**6) - qajp2/(120*dx**6) + qajp3/(720*dx**6)
    b = -qajm1/(48*dx**5) + qajm2/(60*dx**5) - qajm3/(240*dx**5) + qajp1/(48*dx**5) - qajp2/(60*dx**5) + qajp3/(240*dx**5)
    c = 61*qaj/(144*dx**4) - 19*qajm1/(64*dx**4) + 3*qajm2/(32*dx**4) - 5*qajm3/(576*dx**4) - 19*qajp1/(64*dx**4) + 3*qajp2/(32*dx**4) - 5*qajp3/(576*dx**4)
    d = 83*qajm1/(288*dx**3) - 13*qajm2/(72*dx**3) + 7*qajm3/(288*dx**3) - 83*qajp1/(288*dx**3) + 13*qajp2/(72*dx**3) - 7*qajp3/(288*dx**3)
    e = -301*qaj/(192*dx**2) + 229*qajm1/(256*dx**2) - 77*qajm2/(640*dx**2) + 37*qajm3/(3840*dx**2) + 229*qajp1/(256*dx**2) - 77*qajp2/(640*dx**2) + 37*qajp3/(3840*dx**2)
    f = -1891*qajm1/(2304*dx) + 559*qajm2/(2880*dx) - 259*qajm3/(11520*dx) + 1891*qajp1/(2304*dx) - 559*qajp2/(2880*dx) + 259*qajp3/(11520*dx)
    g = 30251*qaj/26880 - 7621*qajm1/107520 + 159*qajm2/17920 - 5*qajm3/7168 - 7621*qajp1/107520 + 159*qajp2/17920 - 5*qajp3/7168
    return a,b,c,d,e,f,g

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

    P3jtojp3a33,P3jtojp3a32,P3jtojp3a31,P3jtojp3a30 = P3jtojp3(qaj,qajp1,qajp2,qajp3,dx)
    P3jm3toja33,P3jm3toja32,P3jm3toja31,P3jm3toja30 =  P3jm3toj(qaj,qajm1,qajm2,qajm3,dx)
    P6jm3tojp3a66,P6jm3tojp3a65,P6jm3tojp3a64,P6jm3tojp3a63,P6jm3tojp3a62,P6jm3tojp3a61,P6jm3tojp3a60  = P6jm3tojp3(qaj,qajm1,qajm2,qajm3,qajp1,qajp2,qajp3,dx)
    
    
    P2jtojp2a22,P2jtojp2a21,P2jtojp2a20 =   P2jtojp2(qaj,qajp1,qajp2,dx)
    P2jm2toja22,P2jm2toja21,P2jm2toja20  =  P2jm2toj(qajm2,qajm1,qaj,dx)
    P4jm2tojp2a44,P4jm2tojp2a43,P4jm2tojp2a42,P4jm2tojp2a41,P4jm2tojp2a40 =  P4jm2tojp2(qaj,qajm1,qajm2,qajp1,qajp2,dx)
    
    P1jtojp1a11,P1jtojp1a10  =  P1jtojp1(qaj,qajp1,dx)
    P1jm1toja11,P1jm1toja10  =  P1jm1toj(qaj,qajm1,dx)
    P2jm1tojp1a22,P2jm1tojp1a21,P2jm1tojp1a20 =  P2jm1tojp1(qajp1,qaj,qajm1,dx)   
    
    Recona3 = minmodL([P3jm3toja33,P6jm3tojp3a63,P3jtojp3a33])
    
    LowRecona2 = minmodL([P2jm2toja22,P4jm2tojp2a42,P2jtojp2a22])
    Recona2 = minmodL([P3jm3toja32-LowRecona2,P6jm3tojp3a62-LowRecona2,P3jtojp3a32-LowRecona2]) + LowRecona2
    
    LowRecona1 = minmodL([P1jm1toja11,P2jm1tojp1a21,P1jtojp1a11])
    Recona1 = minmodL([P3jm3toja31-LowRecona1,P6jm3tojp3a61-LowRecona1,P3jtojp3a31-LowRecona1]) + LowRecona1
    
    Recona0 = qaj - Recona2/3*(dx/2)**2 
    
    # return P3ha66,P3ha65,P3ha64,P3ha63,P3ha62,P3ha61,P3ha60
    # print('Weights',w3,w2)
    # print('H',P3ha66,P3ha65,P3ha64,P3ha63,P3ha62,P3ha61,P3ha60)
    # print('R',a36,a35,a34,na33,na32,na31,na30)
    return  Recona3, Recona2, Recona1, Recona0

def Recon(qA,x,dx,nGcells,eps):
    n = len(x)
    qjphs = zeros(n)
    xjphs = zeros(n)
    qjmhs = zeros(n)
    xjmhs = zeros(n)
    
    for i in range(nGcells,n-nGcells):
        a3,a2,a1,a0 = ReconCoeffs(qA[i-3],qA[i-2],qA[i-1],qA[i],qA[i+1],qA[i+2],qA[i+3],dx)

        #Reconstruction at  x_{j+1/2
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
RightCoeffs = [-a0,a1,a2,a3,a4,a5,a6,a7,a8,a9]

eps = 10.0**(-10)
nGcells = 4
expn = 10
# expn = 1
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