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
        qajm1 = q[i-1]
        qajm2 = q[i-2]
        
        #Start with constant
        P0ja0 = qaj
        
        P0jplot =  0*(xplot - x[i])  +P0ja0

        P1jtojp1a11,P1jtojp1a10  =  P1jtojp1(qaj,qajp1,dx)
        P1jm1toja11,P1jm1toja10  =  P1jm1toj(qaj,qajm1,dx)
        P1jm1tojp1a22,P1jm1tojp1a21,P1jm1tojp1a20  =  P2jm1tojp1(qajp1,qaj,qajm1,dx)
        
        P1Ra11 = minmodL([P1jm1toja11,P1jm1tojp1a21,P1jtojp1a11 ])
        P1Ra10 = qaj
        
        
        P1jtojp1plot =  P1jtojp1a11*(xplot - x[i])  +P1jtojp1a10
        P1jm1tojplot =  P1jm1toja11*(xplot - x[i])  +P1jm1toja10
        LimLinPlot = P1Ra11*(xplot - x[i])  + P1Ra10
        
  
        #Quadratic
        P2jtojp2a22,P2jtojp2a21,P2jtojp2a20=   P2jtojp2(qaj,qajp1,qajp2,dx)
        P2jm2toja22,P2jm2toja21,P2jm2toja20  =  P2jm2toj(qajm2,qajm1,qaj,dx)   
        P4jm2tojp2a44,P4jm2tojp2a43,P4jm2tojp2a42,P4jm2tojp2a41,P4jm2tojp2a40  =  P4jm2tojp2(qaj,qajm1,qajm2,qajp1,qajp2,dx)

        P2Ra22 = minmodL([P2jm2toja22,P4jm2tojp2a42,P2jtojp2a22])   
        P2Ra21 = minmodL([P2jm2toja21,P4jm2tojp2a41,P2jtojp2a21])        
        P2Ra20 = qaj- P2Ra22/3*(dx/2)**2 
        LimP2Plot = P2Ra22*(xplot - x[i])**2 + P2Ra21*(xplot - x[i]) + P2Ra20 
        
        
        P0RLE = abs(P0ja0 - q[i-3])
        P0RRE = abs(P0ja0 - q[i+3])
        
        P1RLE = PolyIntegralError(-7*dx/2,-5*dx/2,q[i-3],dx,(P1Ra10,P1Ra11))
        P1RRE = PolyIntegralError(5*dx/2,7*dx/2,q[i+3],dx,(P1Ra10,P1Ra11))

        P2RLE = PolyIntegralError(-7*dx/2,-5*dx/2,q[i-3],dx,(P2Ra20,P2Ra21,P2Ra22))
        P2RRE = PolyIntegralError(5*dx/2,7*dx/2,q[i+3],dx,(P2Ra20,P2Ra21,P2Ra22))
        
        print(i,xmh,xph,'Recon P0  (j)',P0RLE,P0RRE)
        print(i,xmh,xph,'Recon P1  (j-1,j+1)',P1RLE,P1RRE)
        print(i,xmh,xph,'Recon P2  (j-2,j+2)',P2RLE,P2RRE)

        if (0.5*(P1RLE+P1RRE) > 0.5*(P2RLE+P2RRE)):
            ReconPlot = LimP2Plot
        elif(0.5*(P0RLE+P0RRE) > 0.5*(P1RLE+P1RRE)):
            ReconPlot = LimLinPlot
        else:
            ReconPlot = P0jplot
        
        ReconPlot = LimP2Plot
        
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
            # plot(xplot, P0jplot , '-b',label='j value')   
            # plot(xplot, P1jtojp1plot , '-g',label='j to j+1 value')   
            # plot(xplot, P1jm1tojplot, '-r',label='j-1 to j value')   

            plot(xplot, P0jplot, '-b',label='Recon P1')   
            plot(xplot, LimLinPlot, '-r',label='Recon P1')   
            plot(xplot, LimP2Plot, '-g',label='Recon P2')   
            plot(xplot, ReconPlot, '--r',label='Choose Small')   
        else:
            # plot(xplot, P0jplot , '-b')   
            # plot(xplot, P1jtojp1plot , '-g')   
            # plot(xplot, P1jm1tojplot, '-r')   
            
            plot(xplot, P0jplot, '-b')   
            plot(xplot, LimLinPlot, '-r')   
            plot(xplot, LimP2Plot, '-g')   
            plot(xplot, ReconPlot, '--r')   





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


# def P4jm2tojp2(qaj,qajm1,qajm2,qajp1,qajp2,dx):
#     a =-qaj/(4*dx**4) - qajm1/(6*dx**4) + qajm2/(24*dx**4) - qajp1/(6*dx**4) + qajp2/(24*dx**4)
#     b= qajm1/(6*dx**3) - qajm2/(12*dx**3) - qajp1/(6*dx**3) + qajp2/(12*dx**3)
#     c = -11*qaj/(8*dx**2) + 3*qajm1/(4*dx**2) - qajm2/(16*dx**2) + 3*qajp1/(4*dx**2) - qajp2/(16*dx**2)
#     d = -17*qajm1/(24*dx) + 5*qajm2/(48*dx) + 17*qajp1/(24*dx) - 5*qajp2/(48*dx)
#     e = 1067*qaj/960 - 29*qajm1/480 + 3*qajm2/640 - 29*qajp1/480 + 3*qajp2/640
#     return a,b,c,d,e



# def P3jtojp3(qaj,qajp1,qajp2,qajp3,dx):
#     a=(3*qajp1 - 3*qajp2 + qajp3 - qaj)/(6*dx**3)
#     b=(-5*qajp1 + 4*qajp2 - qajp3 + 2*qaj)/(2*dx**2)
#     c=(69*qajp1 - 33*qajp2 + 7*qajp3 - 43*qaj)/(24*dx)
#     d=5*qajp1/24 - qajp2/6 + qajp3/24 + 11*qaj/12
#     return a,b,c,d

# def P3jm3toj(qaj,qajm1,qajm2,qajm3,dx):
#     a=(-3*qajm1 + 3*qajm2 - qajm3 + qaj)/(6*dx**3)
#     b=(-5*qajm1 + 4*qajm2 - qajm3 + 2*qaj)/(2*dx**2)
#     c=(-69*qajm1 + 33*qajm2 - 7*qajm3 + 43*qaj)/(24*dx)
#     d=5*qajm1/24 - qajm2/6 + qajm3/24 + 11*qaj/12
#     return a,b,c,d

# def P6jm3tojp3(qaj,qajm1,qajm2,qajm3,qajp1,qajp2,qajp3,dx):
#     a =-qaj/(36*dx**6) + qajm1/(48*dx**6) - qajm2/(120*dx**6) + qajm3/(720*dx**6) + qajp1/(48*dx**6) - qajp2/(120*dx**6) + qajp3/(720*dx**6)
#     b = -qajm1/(48*dx**5) + qajm2/(60*dx**5) - qajm3/(240*dx**5) + qajp1/(48*dx**5) - qajp2/(60*dx**5) + qajp3/(240*dx**5)
#     c = 61*qaj/(144*dx**4) - 19*qajm1/(64*dx**4) + 3*qajm2/(32*dx**4) - 5*qajm3/(576*dx**4) - 19*qajp1/(64*dx**4) + 3*qajp2/(32*dx**4) - 5*qajp3/(576*dx**4)
#     d = 83*qajm1/(288*dx**3) - 13*qajm2/(72*dx**3) + 7*qajm3/(288*dx**3) - 83*qajp1/(288*dx**3) + 13*qajp2/(72*dx**3) - 7*qajp3/(288*dx**3)
#     e = -301*qaj/(192*dx**2) + 229*qajm1/(256*dx**2) - 77*qajm2/(640*dx**2) + 37*qajm3/(3840*dx**2) + 229*qajp1/(256*dx**2) - 77*qajp2/(640*dx**2) + 37*qajp3/(3840*dx**2)
#     f = -1891*qajm1/(2304*dx) + 559*qajm2/(2880*dx) - 259*qajm3/(11520*dx) + 1891*qajp1/(2304*dx) - 559*qajp2/(2880*dx) + 259*qajp3/(11520*dx)
#     g = 30251*qaj/26880 - 7621*qajm1/107520 + 159*qajm2/17920 - 5*qajm3/7168 - 7621*qajp1/107520 + 159*qajp2/17920 - 5*qajp3/7168
#     return a,b,c,d,e,f,g


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
q = PB_IC(x,dx)
qh = PB_IC(xh,dx)

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


