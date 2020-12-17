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
        
        P0jplot =  0*(xplot - x[i])  +P0ja0

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
        if TotalSI > 10.0**(-12):
            P1Ra1 = (P1jtojp1SI/TotalSI)*P1jm1toja11 + (P1jm1tojSI/TotalSI)*P1jtojp1a11 
            P1Ra0 = (P1jtojp1SI/TotalSI)*P1jm1toja10 + (P1jm1tojSI/TotalSI)*P1jtojp1a10 
        else:
            P1Ra1 = 0
            P1Ra0 = qaj

   
        LimLinPlot = P1Ra1*(xplot - x[i])  +P1Ra0
         
            
        TotalSI =  (P2jtojp2SI + P2jm1tojp1SI + P2jm2tojSI)*2.0
        if TotalSI > 10.0**(-12):
            P2Ra2 = (((P2jm1tojp1SI + P2jtojp2SI) /TotalSI)*P2jm2toja22 + ((P2jm2tojSI + P2jtojp2SI) /TotalSI)*P2jm1tojp1a22 +  ((P2jm2tojSI + P2jm1tojp1SI) /TotalSI)*P2jtojp2a22)
            P2Ra1 = (((P2jm1tojp1SI + P2jtojp2SI)/TotalSI)*P2jm2toja21 + ((P2jm2tojSI + P2jtojp2SI) /TotalSI)*P2jm1tojp1a21 +  ((P2jm2tojSI + P2jm1tojp1SI) /TotalSI)*P2jtojp2a21)
            P2Ra0 = (((P2jm1tojp1SI + P2jtojp2SI) /TotalSI)*P2jm2toja20 + ((P2jm2tojSI + P2jtojp2SI) /TotalSI)*P2jm1tojp1a20 +  ((P2jm2tojSI + P2jm1tojp1SI)/TotalSI)*P2jtojp2a20 )
        else:
            P2Ra2 = 0
            P2Ra1 = 0
            P2Ra0 = qaj           

        LimP2Plot = P2Ra2*(xplot - x[i])**2 + P2Ra1*(xplot - x[i])   +P2Ra0
        
        
        
        TotalSI =  (Bjm3toj + Bjm2tojp1  +  Bjm1tojp2  +Bjtojp3 )*3.0
        if TotalSI > 10.0**(-12):
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
        
        LimP3Plot = P3Ra3*(xplot - x[i])**2 +P3Ra2*(xplot - x[i])**2 + P3Ra1*(xplot - x[i])   +P3Ra0
        
        # # if (P1Ra1 == 0):
        # #     ReconPlot = P0jplot
        # # elif(P2Ra2 == 0):
        # #     ReconPlot = LimLinPlot
        # # elif(P3Ra3 == 0):
        # #     ReconPlot = LimP2Plot
        # # else:
        ReconPlot = LimP3Plot
            
        P0RE = PolyIntegralError(-7*dx/2,7*dx/2,q[i-3] + q[i-2]+ q[i-1] +q[i] +q[i+1] +q[i+2] + q[i+3],dx,[P0ja0])
        P1RE = PolyIntegralError(-7*dx/2,7*dx/2,q[i-3] + q[i-2]+ q[i-1] +q[i] +q[i+1] +q[i+2] + q[i+3],dx,(P1Ra0,P1Ra1))
        P2RE = PolyIntegralError(-7*dx/2,7*dx/2,q[i-3] + q[i-2]+ q[i-1] +q[i] +q[i+1] +q[i+2] + q[i+3],dx,(P2Ra0,P2Ra1,P2Ra2))
        P3RE = PolyIntegralError(-7*dx/2,7*dx/2,q[i-3] + q[i-2]+ q[i-1] +q[i] +q[i+1] +q[i+2] + q[i+3],dx,(P3Ra0,P3Ra1,P3Ra2,P3Ra3))

        print(i,xmh,xph,'Errors',P0RE,P1RE,P2RE,P3RE)
        
        # ReconPlot = LimP2Plot
        
        # if (max(P0jLE,P0jRE) > max(P1RLE,P1RRE)):
        #     ReconPlot = LimLinPlot
        #     # P0RLE = PolyIntegralError(-7*dx/2,-5*dx/2,q[i-3],dx,[P0ja0])
        #     # P0RRE = PolyIntegralError(5*dx/2,7*dx/2,q[i+3],dx,[P0ja0])
            
        #     # P1RLE = PolyIntegralError(-7*dx/2,-5*dx/2,q[i-3],dx,(P1Ra10,P1Ra11))
        #     # P1RRE = PolyIntegralError(5*dx/2,7*dx/2,q[i+3],dx,(P1Ra10,P1Ra11))
 
        #     # P2RLE = PolyIntegralError(-7*dx/2,-5*dx/2,q[i-3],dx,(P2Ra20,P2Ra21,P2Ra22))
        #     # P2RRE = PolyIntegralError(5*dx/2,7*dx/2,q[i+3],dx,(P2Ra20,P2Ra21,P2Ra22))
            
        #     # if (max(P1RLE,P1RRE) > max(P2RLE,P2RRE)):
        #     #     ReconPlot = LimP2Plot
        #     # elif(max(P0RLE,P0RRE) > max(P1RLE,P1RRE)):
        #     #     ReconPlot = LimLinPlot
        #     # else:
        #     #     ReconPlot = P0jplot
        # else:
        #     ReconPlot = P0jplot
        
        if i == nG:
            plot(xplot, P0jplot, '-b',label='Recon P1')   
            plot(xplot, LimLinPlot, '-r',label='Recon P1')   
            plot(xplot, LimP2Plot, '-g',label='Recon P2')   
            plot(xplot, LimP3Plot, '-y',label='Recon P3')   
            plot(xplot, ReconPlot, '--c',label='Choose Small')   
        else:
            plot(xplot, P0jplot, '-b')   
            plot(xplot, LimLinPlot, '-r')   
            plot(xplot, LimP2Plot, '-g')   
            plot(xplot, LimP3Plot, '-y')  
            plot(xplot, ReconPlot, '--c')   





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
# 
# q = Lin_IC(x)
# qh = Lin_IC(xh)
# 
# q = PB_IC(x,dx)
# qh = PB_IC(xh,dx)

plot(xh,qh,'--k', label='Analytic Value')
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


