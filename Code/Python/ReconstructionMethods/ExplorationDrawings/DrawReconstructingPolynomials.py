from numpy import *
from numpy.linalg import solve,norm
from matplotlib.pyplot import plot,loglog,title,xlabel,ylabel,legend,xlim,ylim



def minmodL(q):
    qmin = min(q)
    qmax = max(q)
    if qmin>0:
        return qmin
    if qmax<0 :
        return qmax
    else:
        return 0


def SP_IC(x,dx):
    
    n = len(x)
    q = zeros(n)
    
    for i in range(n):
        if (x[i] < dx/2 and x[i] > - dx/2):
            q[i] = 1
    
    return q    


def DB_IC(x,dx,ldx):
    
    n = len(x)
    q = ones(n)
    
    for i in range(n):
        if (x[i] < 0.1*ldx):
            q[i] = 2
            
    
    return q  

def EXPPeak_IC(x):
    
    n = len(x)
    q = exp(-abs(x))
    
    return q  


def Lin_IC(x):
    
    n = len(x)
    q = zeros(n)
    for i in range(n):
        if x[i] <= 0:
            q[i] = 20
        else:
            q[i] = 20 - (x[i])
    
    return q


def PB_IC(x,dx):
    
    n = len(x)
    xeh = x[int(n/2)] + 0.5*dx
    q = zeros(n)
    for i in range(n):
        if x[i] <= 0:
            q[i] = 2
        else:
            xmjh = x[i] - dx/2
            xpjh = x[i] + dx/2
            qamjh = 2*xmjh - ((xmjh )**3)/3.0
            qapjh = 2*xpjh - ((xpjh )**3)/3.0
            q[i] = (qapjh - qamjh)/ dx
    
    return q
    
   




def PlotAverages(x,q,dx):
    n = len(x)
    for i in range(n):
        xmh = x[i] - 0.5*dx
        xph = x[i] + 0.5*dx
        plot([xmh,xph], [q[i] , q[i]] , '--k')


          


def PlotLimLinears(x,q,dx,nG,np):
    n = len(x)
    theta  = 1
    for i in range(nG,n-nG):
        xmh = x[i] - 0.5*dx
        xph = x[i] + 0.5*dx
        
        ndx = (xph - xmh)/np
        xplot = arange(xmh,xph + 0.5*ndx ,ndx)
        
        qaj = q[i]
        qajp1 = q[i+1]
        qajm1 = q[i-1]
        P1jtojp1a11,P1jtojp1a10  =  P1jtojp1(qaj,qajp1,dx)
        P1jm1toja11,P1jm1toja10  =  P1jm1toj(qaj,qajm1,dx)
        P2jm1tojp1a22,P2jm1tojp1a21,P2jm1tojp1a20 =  P2jm1tojp1(qajp1,qaj,qajm1,dx)
        
        
        qjphF = P1jtojp1a11*(dx/2) + P1jtojp1a10
        qjmhF = P1jtojp1a11*(-dx/2) + P1jtojp1a10
        
        qjphB = P1jm1toja11*(dx/2) + P1jm1toja10
        qjmhB = P1jm1toja11*(-dx/2) + P1jm1toja10
        
        qjphC = P2jm1tojp1a21*(dx/2) + qaj
        qjmhC = P2jm1tojp1a21*(-dx/2) + qaj
        
        Recona1 = minmodL([P1jm1toja11,P2jm1tojp1a21,P1jtojp1a11])
        Recona0 = qaj
        
        Poly0 = P2jm1tojp1a22*(xplot - x[i])**2 + P2jm1tojp1a21*(xplot - x[i]) +P2jm1tojp1a20
        PolyR = Recona1*(xplot - x[i]) +Recona0
        
        if i == nG:
            plot([xmh,xph], [qjmhF  , qjphF] , '-r',label='Forward')     
            plot([xmh,xph], [qjmhB  , qjphB] , '-b',label='Backward')   
            plot([xmh,xph], [qjmhC  , qjphC] , '-g',label='Centered')   
            plot(xplot, PolyR , '--y',label='Recon')   

        else:
            plot([xmh,xph], [qjmhF  , qjphF] , '-r')     
            plot([xmh,xph], [qjmhB  , qjphB] , '-b')   
            plot([xmh,xph], [qjmhC  , qjphC] , '-g')   
            plot(xplot, PolyR,'--y')  


def PlotLimQuadratics(x,q,dx,nG,np):
    n = len(x)
    theta  = 1
    for i in range(nG,n-nG):
        xmh = x[i] - 0.5*dx
        xph = x[i] + 0.5*dx
        
        ndx = (xph - xmh)/np
        xplot = arange(xmh,xph + 0.5*ndx ,ndx)
        
        qaj = q[i]
        qajp1 = q[i+1]
        qajp2 = q[i+2]
        qajm1 = q[i-1]
        qajm2 = q[i-2]
        
        
        P2jtojp2a22,P2jtojp2a21,P2jtojp2a20 =   P2jtojp2(qaj,qajp1,qajp2,dx)
        P2jm2toja22,P2jm2toja21,P2jm2toja20  =  P2jm2toj(qajm2,qajm1,qaj,dx)
        P4jm2tojp2a44,P4jm2tojp2a43,P4jm2tojp2a42,P4jm2tojp2a41,P4jm2tojp2a40 =  P4jm2tojp2(qaj,qajm1,qajm2,qajp1,qajp2,dx)
        
        P1jtojp1a11,P1jtojp1a10  =  P1jtojp1(qaj,qajp1,dx)
        P1jm1toja11,P1jm1toja10  =  P1jm1toj(qaj,qajm1,dx)
        P2jm1tojp1a22,P2jm1tojp1a21,P2jm1tojp1a20 =  P2jm1tojp1(qajp1,qaj,qajm1,dx)

        
        PolyB = P2jm2toja22*(xplot - x[i])**2 + P2jm2toja21*(xplot - x[i]) +P2jm2toja20
        PolyCH = P4jm2tojp2a42*(xplot - x[i])**2 + P4jm2tojp2a41*(xplot - x[i]) +P4jm2tojp2a40
        PolyCL = P2jm1tojp1a22*(xplot - x[i])**2 + P2jm1tojp1a21*(xplot - x[i]) +P2jm1tojp1a20
        PolyF = P2jtojp2a22*(xplot - x[i])**2 + P2jtojp2a21*(xplot - x[i]) +P2jtojp2a20
        
        Recona2 = minmodL([P2jm2toja22,P4jm2tojp2a42,P2jtojp2a22])
        
        LinearCoeffLim = minmodL([P1jm1toja11,P2jm1tojp1a21,P1jtojp1a11])
        Coeff1At = minmodL([P2jm2toja21 - LinearCoeffLim,P4jm2tojp2a41- LinearCoeffLim,P2jtojp2a21- LinearCoeffLim])  + LinearCoeffLim
        
        Recona1 = Coeff1At
        Recona0 = qaj - Recona2/3*(dx/2)**2 
        
        PolyR = Recona2*(xplot - x[i])**2 + Recona1*(xplot - x[i]) +Recona0
        
        print(i,qaj,Recona2,Recona1,Recona0)
        
        if i == nG:
            plot(xplot, PolyB , '-b',label='Backward')   
            plot(xplot, PolyCH , '-g',label='Centered High')   
            plot(xplot, PolyCL , '-y',label='Centered Low')  
            plot(xplot, PolyF , '-r',label='Forward')   
            plot(xplot, PolyR , '--c',label='Recon Cand')   

        else:
            plot(xplot, PolyB , '-b')   
            plot(xplot, PolyCH , '-g')   
            plot(xplot, PolyCL, '-y')   
            plot(xplot, PolyF , '-r') 
            plot(xplot, PolyR , '--c')   



def PlotLimCubics(x,q,dx,nG,np):
    n = len(x)
    theta  = 1
    for i in range(nG,n-nG):
        xmh = x[i] - 0.5*dx
        xph = x[i] + 0.5*dx
        
        ndx = (xph - xmh)/np
        xplot = arange(xmh,xph + 0.5*ndx ,ndx)
        
        qaj = q[i]
        qajp1 = q[i+1]
        qajp2 = q[i+2]
        qajp3 = q[i+3]
        qajm1 = q[i-1]
        qajm2 = q[i-2]
        qajm3 = q[i-3]
        
        P3jtojp3a33,P3jtojp3a32,P3jtojp3a31,P3jtojp3a30 = P3jtojp3(qaj,qajp1,qajp2,qajp3,dx)
        P3jm3toja33,P3jm3toja32,P3jm3toja31,P3jm3toja30 =  P3jm3toj(qaj,qajm1,qajm2,qajm3,dx)
        P6jm3tojp3a66,P6jm3tojp3a65,P6jm3tojp3a64,P6jm3tojp3a63,P6jm3tojp3a62,P6jm3tojp3a61,P6jm3tojp3a60  = P6jm3tojp3(qaj,qajm1,qajm2,qajm3,qajp1,qajp2,qajp3,dx)
        
        
        P2jtojp2a22,P2jtojp2a21,P2jtojp2a20 =   P2jtojp2(qaj,qajp1,qajp2,dx)
        P2jm2toja22,P2jm2toja21,P2jm2toja20  =  P2jm2toj(qajm2,qajm1,qaj,dx)
        P4jm2tojp2a44,P4jm2tojp2a43,P4jm2tojp2a42,P4jm2tojp2a41,P4jm2tojp2a40 =  P4jm2tojp2(qaj,qajm1,qajm2,qajp1,qajp2,dx)
        
        P1jtojp1a11,P1jtojp1a10  =  P1jtojp1(qaj,qajp1,dx)
        P1jm1toja11,P1jm1toja10  =  P1jm1toj(qaj,qajm1,dx)
        P2jm1tojp1a22,P2jm1tojp1a21,P2jm1tojp1a20 =  P2jm1tojp1(qajp1,qaj,qajm1,dx)

        
        PolyF = P3jtojp3a33*(xplot - x[i])**3 + P3jtojp3a32*(xplot - x[i])**2 + P3jtojp3a31*(xplot - x[i]) + P3jtojp3a30
        PolyB = P3jm3toja33*(xplot - x[i])**3 + P3jm3toja32*(xplot - x[i])**2 + P3jm3toja31*(xplot - x[i]) + P3jm3toja30
        PolyC = P6jm3tojp3a63*(xplot - x[i])**3 + P6jm3tojp3a62*(xplot - x[i])**2 +P6jm3tojp3a61*(xplot - x[i]) + P6jm3tojp3a60
        
        
        Recona3 = minmodL([P3jm3toja33,P6jm3tojp3a63,P3jtojp3a33])
        
        LowRecona2 = minmodL([P2jm2toja22,P4jm2tojp2a42,P2jtojp2a22])
        Recona2 = minmodL([P3jm3toja32-LowRecona2,P6jm3tojp3a62-LowRecona2,P3jtojp3a32-LowRecona2]) + LowRecona2
        
        LowRecona1 = minmodL([P1jm1toja11,P2jm1tojp1a21,P1jtojp1a11])
        Recona1 = minmodL([P3jm3toja31-LowRecona1,P6jm3tojp3a61-LowRecona1,P3jtojp3a31-LowRecona1]) + LowRecona1
        
        Recona0 = qaj - Recona2/3*(dx/2)**2 
        
        PolyR = Recona3*(xplot - x[i])**3 +Recona2*(xplot - x[i])**2 + Recona1*(xplot - x[i]) +Recona0

        
        if i == nG:
            plot(xplot, PolyB , '-b',label='Backward')   
            plot(xplot, PolyC , '-g',label='Centered')   
            plot(xplot, PolyF , '-r',label='Forward')   
            plot(xplot, PolyR , '--c',label='Recon Cand')   

        else:
            plot(xplot, PolyB , '-b')   
            plot(xplot, PolyC , '-g')   
            plot(xplot, PolyF , '-r') 
            plot(xplot, PolyR , '--c')   



#constant
def P0j(q,j):
    return q[j]


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


nG = 5
sx = -2.0
ex = 2.0
dx = 0.5
hdx = dx/20.0
theta = 1.0
x = arange(sx - nG*dx,ex + (nG+1)*dx,dx)
xh = arange(sx- nG*dx,ex + (nG+1)*dx,hdx)


# q =  SP_IC(x,dx) 
# qh =  SP_IC(xh,hdx)

# q =  DB_IC(x,dx,dx) 
# qh =  DB_IC(xh,hdx,dx)

# q =  EXPPeak_IC(x)
# qh =  EXPPeak_IC(xh)

# q = Lin_IC(x)
# qh = Lin_IC(xh)

q = PB_IC(x,dx)
qh = PB_IC(xh,dx)

plot(xh,qh,'-k', label='Analytic Value')
plot(x,q,'.k', label='Average Values')
PlotAverages(x,q,dx)


# PlotLimLinears(x,q,dx,nG,20)
# PlotLimQuadratics(x,q,dx,nG,20)
PlotLimCubics(x,q,dx,nG,20)

legend()
s = "Reconstructions Parabolic Jump"
title(s)
ylabel('q(x)')
xlabel('x')
# xlim((sx,ex))
# ylim((min(q[nG :-nG])*0.99,max(q[nG :-nG])*1.01))

#Look at different Linears


