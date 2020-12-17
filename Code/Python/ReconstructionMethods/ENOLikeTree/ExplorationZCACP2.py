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
    q = exp(-(x**2))
    q =4 -(x**3)
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


def PlotBuild1(x,q,dx,nG,np):
    n = len(x)
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
        
        if (abs(P1jtojp1SI -P1jm1tojSI) < 10.0**(-14)):
            I1 = 0
        #To left is more smooth
        elif P1jtojp1SI > P1jm1tojSI:
            I1 = -1
            R2a1 = P1jm1toja11
            R2a0 = P1jm1toja10

        #To right is more smooth
        else:
            I1 = 1
            R2a1 = P1jtojp1a11 
            R2a0 = P1jtojp1a10 
            

        if I1 ==0:
            if(P1jtojp1a11 == P1jm1toja11 and P1jtojp1a10 == P1jm1toja10 ):
                R2a1 = P1jtojp1a11
                R2a0 = P1jtojp1a10
            else:
                R2a1 = 0
                R2a0 = qaj            
        R2a2 = 0
        LimP2R = R2a2*(xplot - x[i])**2 + R2a1*(xplot - x[i])   +R2a0

        
        RIxjph = R2a2/3*(dx/2)**3 + R2a1/2*(dx/2)**2 + R2a0*(dx/2)
        RIxjmh = R2a2/3*(-dx/2)**3 + R2a1/2*(-dx/2)**2 + R2a0*(-dx/2)
        print(i,xmh,xph,'CA', (RIxjph -  RIxjmh)/dx,qaj,(RIxjph -  RIxjmh)/dx-qaj )
        
        
        P0jplot =  0*(xplot - x[i])  +P0ja0      
        LimP1Plotjm1j = P1jm1toja11*(xplot - x[i])   +P1jm1toja10
        LimP1Plotjjp1 = P1jtojp1a11*(xplot - x[i])   +P1jtojp1a10
        
        if i == nG:            
            plot(xplot, LimP2R, '-b',label='Recon P2 j-2,j+2')
        else:
            plot(xplot, LimP2R, '-b')
            


def PlotBuild2(x,q,dx,nG,np):
    n = len(x)
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
        
        if (abs(P1jtojp1SI -P1jm1tojSI) < 10.0**(-14)):
            I1 = 0
        #To left is more smooth
        elif P1jtojp1SI > P1jm1tojSI:
            I1 = -1
            Ra1 = P1jm1toja11
            Ra0 = P1jm1toja10

        #To right is more smooth
        else:
            I1 = 1
            Ra1 = P1jtojp1a11 
            Ra0 = P1jtojp1a10 
            

        # Quadratic Terms
        P2jm2toja22,P2jm2toja21,P2jm2toja20,P2jm2tojSI = P2jm2toj(qajm2,qajm1,qaj,dx)   
        P2jm1tojp1a22,P2jm1tojp1a21,P2jm1tojp1a20,P2jm1tojp1SI  = P2jm1tojp1(qajp1,qaj,qajm1,dx)
        P2jtojp2a22,P2jtojp2a21,P2jtojp2a20,P2jtojp2SI =   P2jtojp2(qaj,qajp1,qajp2,dx) 
        
        #Right
        if I1 == 1:
            if(P2jm1tojp1SI > P2jtojp2SI ):
                Ra2 = P2jtojp2a22
                Ra1 = Ra1 - dx*Ra2
                Ra0 = Ra0 - dx**2/12*Ra2
            else:
                I1 = 0
        #Left
        elif I1 == -1:
            if(P2jm1tojp1SI > P2jm2tojSI ):
                Ra2 = P2jm2toja22
                Ra1 = Ra1 + dx*Ra2
                Ra0 = Ra0 - dx**2/12*Ra2
            else:
                I1 = 0
        #undecided
        if (I1 == 0):
     
            if( P2jm1tojp1SI < P2jm2tojSI or  P2jm1tojp1SI < P2jtojp2SI  ):
                Ra2 = P2jm1tojp1a22
                Ra1 = P2jm1tojp1a21
                Ra0 = P2jm1tojp1a20
            elif( abs(P1jtojp1SI- P1jm1tojSI) < 10.0**(-14)  and P1jtojp1a11 == P1jm1toja11 and P1jtojp1a10 == P1jm1toja10 ):
                Ra2 = 0
                Ra1 = P1jtojp1a11
                Ra0 = P1jtojp1a10
            else:
                Ra2 =0
                Ra1 = 0
                Ra0 = qaj  

        print (i,xmh,xph,'SI 1',P1jtojp1SI,P1jm1tojSI, abs(P1jtojp1SI-P1jm1tojSI)  )      
        print (i,xmh,xph,'SI 2',P2jm2tojSI,P2jm1tojp1SI, P2jtojp2SI,abs(P2jm2tojSI-P2jm1tojp1SI),  abs(P2jm2tojSI -P2jtojp2SI ), abs(P2jm1tojp1SI- P2jtojp2SI) )    
        
        LimP2R = Ra2*(xplot - x[i])**2 + Ra1*(xplot - x[i])   +Ra0

        
        RIxjph = Ra2/3*(dx/2)**3 + Ra1/2*(dx/2)**2 + Ra0*(dx/2)
        RIxjmh = Ra2/3*(-dx/2)**3 + Ra1/2*(-dx/2)**2 + Ra0*(-dx/2)
        print(i,xmh,xph,'CA', (RIxjph -  RIxjmh)/dx,qaj,(RIxjph -  RIxjmh)/dx-qaj )
        
        
        P0jplot =  0*(xplot - x[i])  +P0ja0      
        LimP1Plotjm1j = P1jm1toja11*(xplot - x[i])   +P1jm1toja10
        LimP1Plotjjp1 = P1jtojp1a11*(xplot - x[i])   +P1jtojp1a10
        
        if i == nG:          
            plot(xplot, LimP2R, '-g',label='Recon P2 j-2,j+2')
        else:
            plot(xplot, LimP2R, '-g')






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

q =  EXPPeak_IC(x)
qh =  EXPPeak_IC(xh)

# q =  EXPPeakS_IC(x)
# qh =  EXPPeakS_IC(xh)

# q = Lin_IC(x)
# qh = Lin_IC(xh)

# q = PB_IC(x,dx)
# qh = PB_IC(xh,dx)

plot(xh,qh,'--k', label='Analytic Value')
plot(x,q,'.k', label='Average Values')
PlotAverages(x,q,dx)


# PlotLimLinears(x,q,dx,nG,20)
# PlotLimQuadratics(x,q,dx,nG,20)
# PlotLimCubics(x,q,dx,nG,20)

PlotBuild1(x,q,dx,nG,20)
PlotBuild2(x,q,dx,nG,20)

legend()
s = "Reconstructions Parabolic Jump"
title(s)
ylabel('q(x)')
xlabel('x')
# xlim((sx,ex))
# ylim((min(q[nG :-nG])*0.99,max(q[nG :-nG])*1.01))

#Look at different Linears


