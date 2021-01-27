from numpy import *
from numpy.linalg import solve,norm
from matplotlib.pyplot import plot,loglog,title,xlabel,ylabel,legend,xlim,ylim

#get list of polynomials (list of list of coefficients [[a0,a1,...],[a0,a1]]) - going from left to right
def ModSign(q):
    eps = 10.0**(-12)
    if abs(q) < eps:
        return 0
    else:
        return sign(q)

def LimP1(Pjm1j,Pjjp1):
    Qq = [Pjm1j,Pjjp1]
    if (CompareCoeffsSign(Qq)):
        return 0.5*(Pjm1j[0] + Pjjp1[0])
    else:
        return 0
    

def CompareCoeffsSign(Qq):
    m = len(Qq[0])
    AllCoeffsSameSign = True
    for i in range(m):
        AllCoeffsSameSign = AllCoeffsSameSign and PLCSign(Qq,i)
    return AllCoeffsSameSign

def PLCSign(Qq,j):
    n = len(Qq)
    Sign = ModSign(Qq[0][j])
    SameSign = True
    for i in range(1,n):
        SameSign = SameSign and (Sign == ModSign(Qq[i][j]))    
    return SameSign        
    

def minmodL(q):
    qmin = min(q)
    qmax = max(q)
    if qmin>0:
        return qmin
    if qmax<0 :
        return qmax
    else:
        return 0

def SAL(q):
    qmin = min(q)
    qmax = max(q)
    if qmin>0:
        return 1
    if qmax<0 :
        return 1
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

def EXPPeakS_IC(x):
    
    n = len(x)
    # q = exp(-(x**2))
    q = 4 -(x**4)
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


def PolyIntegralError(LE,RE,AV,dx,Coeffs):
    LEInt = 0
    REInt=0
    n = len(Coeffs)
    
    for i in range(n):
        LEInt = LEInt + Coeffs[i]/(i+1)*LE**(i+1)
        REInt = REInt + Coeffs[i]/(i+1)*RE**(i+1)
    CA = (REInt - LEInt)/dx
    return abs(CA - AV)

def PlotBuildUpToQuad(x,q,dx,nG,np):
    n = len(x)
    errtol = 10.0**(-14)
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
        
        #Start with constant
        P0ja0 = qaj

        #Linear
        P1jtojp1a11,P1jtojp1a10,P1jtojp1SI  =  P1jtojp1(qaj,qajp1,dx)
        P1jm1toja11,P1jm1toja10,P1jm1tojSI  =  P1jm1toj(qaj,qajm1,dx)

            
            
        # Quadratic Terms
        P2jm2toja22,P2jm2toja21,P2jm2toja20,P2jm2tojSI = P2jm2toj(qajm2,qajm1,qaj,dx)   
        P2jm1tojp1a22,P2jm1tojp1a21,P2jm1tojp1a20,P2jm1tojp1SI  = P2jm1tojp1(qajp1,qaj,qajm1,dx)
        P2jtojp2a22,P2jtojp2a21,P2jtojp2a20,P2jtojp2SI =   P2jtojp2(qaj,qajp1,qajp2,dx) 
        
        P4jm2tojp2a44,P4jm2tojp2a43,P4jm2tojp2a42,P4jm2tojp2a41,P4jm2tojp2a40,P4jm2tojp22SI  = P4jm2tojp2(qajm2,qajm1,qaj,qajp1,qajp2,dx)


        #Cubics Terms
        P3jm3toja33,P3jm3toja32,P3jm3toja31,P3jm3toja30,P3jm3tojSI = P3jm3toj(qaj,qajm1,qajm2,qajm3,dx)
        P3jm2tojp1a33,P3jm2tojp1a32,P3jm2tojp1a31,P3jm2tojp1a30,P3jm2tojp1SI = P3jm2tojp1(qaj,qajm1,qajm2,qajp1,dx)
        P3jm1tojp2a33,P3jm1tojp2a32,P3jm1tojp2a31,P3jm1tojp2a30,P3jm1tojp2SI = P3jm1tojp2(qaj,qajm1,qajp1,qajp2,dx)
        P3jtojp3a33,P3jtojp3a32,P3jtojp3a31,P3jtojp3a30,P3jtojp3SI = P3jtojp3(qaj,qajp1,qajp2,qajp3,dx)
        
        P6jm3tojp3a66,P6jm3tojp3a65,P6jm3tojp3a64,P6jm3tojp3a63,P6jm3tojp3a62,P6jm3tojp3a61,P6jm3tojp3a60,P6jm3tojp3SI =P6jm3tojp3(qajm3,qajm2,qajm1,qaj,qajp1,qajp2,qajp3,dx)
        
        P0jP  = 0*(xplot - x[i]) + P0ja0 
        
        P1jtojp1P =P1jtojp1a11*(xplot - x[i]) + P1jtojp1a10
        P1jm1tojP =P1jm1toja11*(xplot - x[i]) + P1jm1toja10
        
        CAP2 = P2jm1tojp1a21*(xplot - x[i])   +qaj
        
        P2jtojp2P = P2jtojp2a22*(xplot - x[i])**2 + P2jtojp2a21*(xplot - x[i])   +P2jtojp2a20
        P2jm1tojp1P = P2jm1tojp1a22*(xplot - x[i])**2 + P2jm1tojp1a21*(xplot - x[i])   +P2jm1tojp1a20
        P2jm2tojP = P2jm2toja22*(xplot - x[i])**2 + P2jm2toja21*(xplot - x[i])   +P2jm2toja20
        
        P4a2 = P4jm2tojp2a42
        P4a1 = P4jm2tojp2a41
        P4a0 = P4jm2tojp2a40
        CAP4 = P4a2*(xplot - x[i])**2 + P4a1*(xplot - x[i])   +P4a0
        

        P3jtojp3P = P3jtojp3a33*(xplot - x[i])**3 + P3jtojp3a32*(xplot - x[i])**2 + P3jtojp3a31*(xplot - x[i])   +P3jtojp3a30
        P3jm1tojp2P = P3jm1tojp2a33*(xplot - x[i])**3 + P3jm1tojp2a32*(xplot - x[i])**2 + P3jm1tojp2a31*(xplot - x[i])   +P3jm1tojp2a30
        P3jm2tojp1P = P3jm2tojp1a33*(xplot - x[i])**3 + P3jm2tojp1a32*(xplot - x[i])**2 + P3jm2tojp1a31*(xplot - x[i])   +P3jm2tojp1a30
        P3jm3tojP = P3jm3toja33*(xplot - x[i])**3 + P3jm3toja32*(xplot - x[i])**2 + P3jm3toja31*(xplot - x[i])   +P3jm3toja30

        P6a3 = P6jm3tojp3a63
        P6a2 = P6jm3tojp3a62
        P6a1 = P6jm3tojp3a61
        P6a0 = P6jm3tojp3a60
        CAP6 = P6a3 *(xplot - x[i])**3 + P6a2 *(xplot - x[i])**2  + P6a1*(xplot - x[i])   +P6a0
        
        print(i,xmh,xph, 'Cells Polynomials')
        print('First Order')
        print('j-1,j   %2.12f %2.12f' %(P1jm1toja10,P1jm1toja11))
        print('j,j+1   %2.12f %2.12f' %(P1jtojp1a10,P1jtojp1a11))
        print()
        print('j-1,j+1 %2.12f %2.12f'%(qaj,P2jm1tojp1a21))
        
        SSIP1 = CompareCoeffsSign([[P1jm1toja11],[P1jtojp1a11]])
        
        if  SSIP1:
            SSIP1V = 0.99*qaj
        else:
            SSIP1V = 1.01*qaj
        
        SSIP2 = CompareCoeffsSign([[P2jm2toja22,P2jm2toja21],[P2jm1tojp1a22,P2jm1tojp1a21],[P2jtojp2a22,P2jtojp2a21]])
        SSIP2rjm2 = CompareCoeffsSign([[P2jm1tojp1a22,P2jm1tojp1a21],[P2jtojp2a22,P2jtojp2a21]])
        SSIP2rjp2 = CompareCoeffsSign([[P2jm2toja22,P2jm2toja21],[P2jm1tojp1a22,P2jm1tojp1a21]])
        if  SSIP2:
            SSIP2V = 0.98*qaj
        else:
            SSIP2V = 1.02*qaj
            
        if  SSIP2rjm2:
            SSIP2rjm2V = 0.97*qaj
        else:
            SSIP2rjm2V = 1.03*qaj
            
        if  SSIP2rjp2:
            SSIP2rjp2V = 0.96*qaj
        else:
            SSIP2rjp2V = 1.04*qaj
        print('Second Order')
        print('j-2,j   %2.12f %2.12f %2.12f'%(P2jm2toja20,P2jm2toja21,P2jm2toja22))
        print('j-1,j+1 %2.12f %2.12f %2.12f'%(P2jm2toja20,P2jm1tojp1a21,P2jm1tojp1a22))
        print('j,j+2   %2.12f %2.12f %2.12f'%(P2jtojp2a20,P2jtojp2a21,P2jtojp2a22))
        print()
        print('j-2,j+2 %2.12f %2.12f %2.12f'%(qaj- P4a2/3*(dx/2)**2   ,P4a1,P4a2))
        
        # # print('Third Order')
        # # print('j-3,j   %2.12f %2.12f %2.12f %2.12f'%(P3jm3toja30,P3jm3toja31,P3jm3toja32,P3jm3toja33))
        # # print('j-2,j+1 %2.12f %2.12f %2.12f %2.12f'%(P3jm2tojp1a30,P3jm2tojp1a31,P3jm2tojp1a32,P3jm2tojp1a33))
        # # print('j-1,j+2 %2.12f %2.12f %2.12f %2.12f'%(P3jm1tojp2a30,P3jm1tojp2a31,P3jm1tojp2a32,P3jm1tojp2a33))
        # # print('j,j+3   %2.12f %2.12f %2.12f %2.12f'%(P3jtojp3a30,P3jtojp3a31,P3jtojp3a32,P3jtojp3a33))
        # # print()
        # # print('j-3,j+3 %2.12f %2.12f %2.12f %2.12f'%(qaj - P6a2/3*(dx/2)**2 ,P6a1,P6a2,P6a3))
        # # print('')


        lw =2
        if i == nG:
            # plot(xplot, P0jP, '-',color ='green',label='Recon P0', linewidth= lw)   
            
            # plot(xplot, P1jm1tojP, '-',color = 'blue',label='P1 j,j-1', linewidth= lw)
            # plot(xplot, P1jtojp1P, '-',color ='red',label='P1 j,j+1', linewidth= lw)
            
            # plot(xplot,  CAP2, '--',color ='lime',label='CP1 j-1,j+1', linewidth=lw)
            
            plot(x[i],SSIP1V,'*r')
            plot(x[i],SSIP2V,'*b')
            plot(x[i],SSIP2rjm2V,'*g')
            plot(x[i],SSIP2rjp2V,'*y')
            
            
            # plot(xplot, P2jm2tojP, '-',color = 'mediumblue',label='P2 j-2,j', linewidth=lw)
            # plot(xplot, P2jm1tojp1P, '-',color ='lime',label='P2 j-1,j+1', linewidth=lw)        
            # plot(xplot, P2jtojp2P, '-',color ='coral',label='P2 j,j+2', linewidth=lw)
            
            # plot(xplot,  CAP4, '--',color ='springgreen',label='CP2 j-2,j+2', linewidth=lw)

            # plot(xplot,P3jm3tojP ,'-',color = 'darkblue',label='P3 j-3,j', linewidth=lw)
            # plot(xplot, P3jm2tojp1P, '-',color ='seagreen',label='P3 j-2,j+1', linewidth=lw)
            # plot(xplot, P3jm1tojp2P, '-',color ='gold',label='P3 j-1,j+2', linewidth=lw)    
            # plot(xplot, P3jtojp3P, '-',color ='brown',label='P3 j,j+3', linewidth=lw) 
            
            # plot(xplot,  CAP6, '--',color ='darkseagreen',label='CP6 j-3,j+3', linewidth=lw)
            
            
        else:
            # plot(xplot, P0jP, '-',color ='green', linewidth= lw)   
            
            # plot(xplot, P1jm1tojP, '-',color = 'blue', linewidth= lw)
            # plot(xplot, P1jtojp1P, '-',color ='red', linewidth= lw)
            
            # plot(xplot,  CAP2, '--',color ='lime', linewidth=lw)
            
            plot(x[i],SSIP1V,'*r')
            plot(x[i],SSIP2V,'*b')
            plot(x[i],SSIP2rjm2V,'*g')
            plot(x[i],SSIP2rjp2V,'*y')           
            
            # plot(xplot, P2jm2tojP, '-',color = 'mediumblue', linewidth=lw)
            # plot(xplot, P2jm1tojp1P, '-',color ='lime', linewidth=lw)        
            # plot(xplot, P2jtojp2P, '-',color ='coral', linewidth=lw)
            
            # plot(xplot,  CAP4, '--',color ='springgreen', linewidth=lw)

            # plot(xplot,P3jm3tojP ,'-',color = 'darkblue', linewidth=lw)
            # plot(xplot, P3jm2tojp1P, '-',color ='lightseagreen', linewidth=lw)
            # plot(xplot, P3jm1tojp2P, '-',color ='gold', linewidth=lw)    
            # plot(xplot, P3jtojp3P, '-',color ='brown', linewidth=lw) 
            
            # plot(xplot,  CAP6, '--',color ='darkseagreen', linewidth=lw)
            


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

nG = 5
sx = -2.0
ex = 2.0
dx = 0.25
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

# q =  EXPPeakS_IC(x)
# qh =  EXPPeakS_IC(xh)

# q = Lin_IC(x)
# qh = Lin_IC(xh)

q = PB_IC(x,dx)
qh = PB_IC(xh,hdx)

# plot(xh,qh,'--k', label='Analytic Value')
plot(x,q,'.k', label='Average Values')
# PlotAverages(x,q,dx)


# PlotLimLinears(x,q,dx,nG,20)
# PlotLimQuadratics(x,q,dx,nG,20)
# PlotLimCubics(x,q,dx,nG,20)

PlotBuildUpToQuad(x,q,dx,nG,20)

legend()
s = "Reconstructions Parabolic Jump"
title(s)
ylabel('q(x)')
xlabel('x')
# xlim((sx,ex))
# ylim((min(q[nG :-nG])*0.99,max(q[nG :-nG])*1.01))

#Look at different Linears


