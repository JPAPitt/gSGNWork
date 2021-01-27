
#Assume we have the expressions at all the necessary points then solved the matrix equation

from numpy import *
from numpy.linalg import solve,norm
from matplotlib.pyplot import plot,loglog


def ModSign(q):
    eps = 10.0**(-12)
    if abs(q) < eps:
        return 0
    else:
        return sign(q)


def LimP2(Pjm2j,Pjm1jp1,Pjjp2,Pjm2jp2,Pjm1j,Pjjp1):
    
    QqFull = [Pjm2j,Pjm1jp1,Pjjp2]
    Qqnojm2 = [Pjm1jp1,Pjjp2]
    Qqnojp2 = [Pjm2j,Pjm1jp1]
    QqFull1 =[Pjm2j,Pjm1jp1,Pjjp2,Pjm2jp2]
    
    Ra2 = 0
    Ra1 = 0
    # PLCSign(Qq,1)
    if (PLCSign(QqFull,1)):
        Ra2,ind1 = minmodL([Pjm2j[1],Pjm1jp1[1],Pjjp2[1],Pjm2jp2[1]] )
        Ra1 = QqFull1[ind1][0]
    else:        
        if (PLCSign(Qqnojm2,1)):
            Ra2,ind1  = minmodL([Pjm1jp1[1],Pjjp2[1]] )
            Ra1 = [Pjm1jp1[0],Pjjp2[0]][ind1]
        elif (PLCSign(Qqnojp2,1)):
            Ra2,ind1  = minmodL([Pjm2j[1],Pjm1jp1[1]] )
            Ra1 = [Pjm2j[0],Pjm1jp1[0]][ind1]
        else:
            Ra1,ia1 = minmodL([Pjm1j[0],Pjjp1[0],0.5*(Pjm1j[0]+Pjjp1[0])])
    return Ra1,Ra2


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
    
    min_index = q.index(qmin)
    max_index = q.index(qmax)
    if qmin>0:
        return qmin,min_index
    if qmax<0 :
        return qmax,max_index
    else:
        return 0,-1

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

def ReconjBE(qajm3,qajm2,qajm1,qaj,qajp1,qajp2,qajp3,dx,eps,j):
    
    #Linear
    #linear
    P1jtojp1a11,P1jtojp1a10,P1jtojp1SI  =  P1jtojp1(qaj,qajp1,dx)
    P1jm1toja11,P1jm1toja10,P1jm1tojSI  =  P1jm1toj(qaj,qajm1,dx)
    
    
    P2jm2toja22,P2jm2toja21,P2jm2toja20,P2jm2tojSI = P2jm2toj(qajm2,qajm1,qaj,dx)   
    P2jm1tojp1a22,P2jm1tojp1a21,P2jm1tojp1a20,P2jm1tojp1SI  = P2jm1tojp1(qajp1,qaj,qajm1,dx)
    P2jtojp2a22,P2jtojp2a21,P2jtojp2a20,P2jtojp2SI =   P2jtojp2(qaj,qajp1,qajp2,dx) 
    
    P4jm2tojp2a44,P4jm2tojp2a43,P4jm2tojp2a42,P4jm2tojp2a41,P4jm2tojp2a40,P4jm2tojp22SI  = P4jm2tojp2(qajm2,qajm1,qaj,qajp1,qajp2,dx)
    
    Ra1,Ra2 = LimP2([P2jm2toja21,P2jm2toja22],[P2jm1tojp1a21,P2jm1tojp1a22],[P2jtojp2a21,P2jtojp2a22],[P4jm2tojp2a41,P4jm2tojp2a42],[P1jm1toja11],[P1jtojp1a11])
    Ra0 = qaj - Ra2/3*(dx/2)**2
    # Ra1 = LimP1([P1jm1toja11],[P1jtojp1a11])
    # Ra0 = qaj

    return Ra2*(-dx/2)**2 + Ra1*(-dx/2) + Ra0, Ra2*(dx/2)**2 + Ra1*(dx/2) + Ra0

# def Reconjph(qajm3,qajm2,qajm1,qaj,qajp1,qajp2,qajp3,dx,eps):

#     pjm3toja,pjm3tojb,pjm3tojc,pjm3tojd,Bjm3toj = P3jm3toj(qaj,qajm1,qajm2,qajm3,dx)
#     pjm2tojp1a,pjm2tojp1b,pjm2tojp1c,pjm2tojp1d,Bjm2tojp1 = P3jm2tojp1(qaj,qajm1,qajm2,qajp1,dx)
#     pjm1tojp2a,pjm1tojp2b,pjm1tojp2c,pjm1tojp2d,Bjm1tojp2 = P3jm1tojp2(qaj,qajm1,qajp1,qajp2,dx)
#     pjtojp3a,pjtojp3b,pjtojp3c,pjtojp3d,Bjtojp3 = P3jtojp3(qaj,qajp1,qajp2,qajp3,dx)
    
            
#     iw1 = ((1.0/35.0) / (eps +Bjm3toj )**2)
#     iw2 = ((12.0/35.0) / (eps +Bjm2tojp1 )**2)
#     iw3 = ((18.0/35.0)  / (eps +Bjm1tojp2 )**2)
#     iw4 = ((4.0/35.0)  / (eps +Bjtojp3 )**2)

    
#     w1 = iw1 / (iw1 + iw2 + iw3 + iw4 )
#     w2 = iw2 / (iw1 + iw2 + iw3 + iw4 )
#     w3 = iw3 / (iw1 + iw2 + iw3 + iw4 )
#     w4 = iw4 / (iw1 + iw2 + iw3 + iw4 )
    
#     qa = w1*pjm3toja + w2*pjm2tojp1a + w3*pjm1tojp2a + w4*pjtojp3a
#     qb = w1*pjm3tojb + w2*pjm2tojp1b + w3*pjm1tojp2b + w4*pjtojp3b
#     qc = w1*pjm3tojc + w2*pjm2tojp1c + w3*pjm1tojp2c + w4*pjtojp3c
#     qd = w1*pjm3tojd + w2*pjm2tojp1d + w3*pjm1tojp2d + w4*pjtojp3d
        

#     return qa*(dx/2)**3 + qb*(dx/2)**2 + qc*(dx/2) + qd


# def Reconjmh(qajm3,qajm2,qajm1,qaj,qajp1,qajp2,qajp3,dx,eps):

#     pjm3toja,pjm3tojb,pjm3tojc,pjm3tojd,Bjm3toj = P3jm3toj(qaj,qajm1,qajm2,qajm3,dx)
#     pjm2tojp1a,pjm2tojp1b,pjm2tojp1c,pjm2tojp1d,Bjm2tojp1 = P3jm2tojp1(qaj,qajm1,qajm2,qajp1,dx)
#     pjm1tojp2a,pjm1tojp2b,pjm1tojp2c,pjm1tojp2d,Bjm1tojp2 = P3jm1tojp2(qaj,qajm1,qajp1,qajp2,dx)
#     pjtojp3a,pjtojp3b,pjtojp3c,pjtojp3d,Bjtojp3 = P3jtojp3(qaj,qajp1,qajp2,qajp3,dx)
    
#     iw1 = ((4.0/35.0) / (eps +Bjm3toj )**2)
#     iw2 = ((18.0/35.0) / (eps +Bjm2tojp1 )**2)
#     iw3 = ((12.0/35.0)  / (eps +Bjm1tojp2 )**2)
#     iw4 = ((1.0/35.0)  / (eps +Bjtojp3 )**2)
    
    
#     w1 = iw1 / (iw1 + iw2 + iw3 + iw4 )
#     w2 = iw2 / (iw1 + iw2 + iw3 + iw4 )
#     w3 = iw3 / (iw1 + iw2 + iw3 + iw4 )
#     w4 = iw4 / (iw1 + iw2 + iw3 + iw4 )
    
#     qa = w1*pjm3toja + w2*pjm2tojp1a + w3*pjm1tojp2a + w4*pjtojp3a
#     qb = w1*pjm3tojb + w2*pjm2tojp1b + w3*pjm1tojp2b + w4*pjtojp3b
#     qc = w1*pjm3tojc + w2*pjm2tojp1c + w3*pjm1tojp2c + w4*pjtojp3c
#     qd = w1*pjm3tojd + w2*pjm2tojp1d + w3*pjm1tojp2d + w4*pjtojp3d
        

#     return qa*(-dx/2)**3 + qb*(-dx/2)**2 + qc*(-dx/2) + qd

def FVMSWWE(ha,uha,nGcells,g,dx,dt,eps):
    
    n = len(ha)
    hap = zeros(n)
    uhap = zeros(n)
    
    j = nGcells-1
    
    hil,hir = ReconjBE(ha[j-3],ha[j-2],ha[j-1],ha[j],ha[j+1],ha[j+2],ha[j+3],dx,eps,j)
    hip1l,hip1r = ReconjBE(ha[j-2],ha[j-1],ha[j],ha[j+1],ha[j+2],ha[j+3],ha[j+4],dx,eps,j)
    
    uhil,uhir = ReconjBE(uha[j-3],uha[j-2],uha[j-1],uha[j],uha[j+1],uha[j+2],uha[j+3],dx,eps,j)
    uhip1l,uhip1r = ReconjBE(uha[j-2],uha[j-1],uha[j],uha[j+1],uha[j+2],uha[j+3],uha[j+4],dx,eps,j)
    
    uir = uhir/ hir
    uip1l = uhip1l/ hip1l
    
    sl = min(0,uir - sqrt(g*hir), uip1l  - sqrt(g*hip1l))
    sr = max(0,uir + sqrt(g*hir), uip1l + sqrt(g*hip1l))
    
    
    felh = uhir
    ferh = uhip1l
    
    feluh = uir*uhir + g*hir*hir/2.0
    feruh = uip1l*uhip1l + g*hip1l*hip1l/2.0
    
    if sr == sl:
        isrmsl = 0
    else:
        isrmsl = 1.0/(sr - sl)
    
    foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hip1l - hir))
    fouh = isrmsl*(sr*feluh - sl*feruh + sl*sr*(uhip1l - uhir))
    

    fih = foh
    fiuh = fouh
    for j in range(nGcells,n-  nGcells):

        hil,hir = ReconjBE(ha[j-3],ha[j-2],ha[j-1],ha[j],ha[j+1],ha[j+2],ha[j+3],dx,eps,j)
        hip1l,hip1r = ReconjBE(ha[j-2],ha[j-1],ha[j],ha[j+1],ha[j+2],ha[j+3],ha[j+4],dx,eps,j)
        
        uhil,uhir = ReconjBE(uha[j-3],uha[j-2],uha[j-1],uha[j],uha[j+1],uha[j+2],uha[j+3],dx,eps,j)
        uhip1l,uhip1r = ReconjBE(uha[j-2],uha[j-1],uha[j],uha[j+1],uha[j+2],uha[j+3],uha[j+4],dx,eps,j)
        
        uir = uhir/ hir
        uip1l = uhip1l/ hip1l
        
        sl = min(0,uir - sqrt(g*hir), uip1l  - sqrt(g*hip1l))
        sr = max(0,uir + sqrt(g*hir), uip1l + sqrt(g*hip1l))
        
        
        felh = uhir
        ferh = uhip1l
        
        feluh = uir*uhir + g*hir*hir/2.0
        feruh = uip1l*uhip1l + g*hip1l*hip1l/2.0
        
        if sr == sl:
            isrmsl = 0
        else:
            isrmsl = 1.0/(sr - sl)
        
        foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hip1l - hir))
        fouh = isrmsl*(sr*feluh - sl*feruh + sl*sr*(uhip1l - uhir))
        
        hap[j] =ha[j] - dt*(foh - fih)/dx
        uhap[j] =uha[j] - dt*(fouh - fiuh)/dx
        
        fih = foh
        fiuh = fouh
    
    for j in range(nGcells):
        hap[j] = ha[j]
        hap[n-1 - j]  = ha[n-1 - j]
        uhap[j] = uha[j]
        uhap[n-1 - j]  = uha[n-1 - j]
        
    return hap,uhap

def EulerStep(hA,uhA,g,dx,nGcells,dt,eps):
    
    hap,uhap = FVMSWWE(hA,uhA,nGcells,g,dx,dt,eps)
    return hap,uhap

def RKStep(hA,uhA,g,dx,nGcells,dt,eps):
    
    hAp,uhAp = EulerStep(hA,uhA,g,dx,nGcells,dt,eps)
    hApp,uhApp = EulerStep(hAp,uhAp,g,dx,nGcells,dt,eps)
    
    hA = 0.5*(hApp + hA)
    uhA = 0.5*(uhApp + uhA)
    return hA ,uhA
    

def Solver(hA,uhA,st,et,g,dx,dt,nGcells,eps):
    
    ct = st
    while ct < st + et:
        
        hA,uhA = RKStep(hA,uhA,g,dx,nGcells,dt,eps)
        # print(hA[len(hA) //2 -1],hA[len(hA) //2 ], hA[len(hA) //2 +1])
        ct = ct + dt
        print(ct)

    return hA,uhA


def initialhDB(x,dx,h0,h1):
    n = len(x)
    hA = zeros(n)
    uhA = zeros(n)

    for i in range(n):
        
        if x[i] < 0:
            hA[i] = h0
        else:
            hA[i] = h1
     
    return hA,uhA

eps = 10.0**(-12)
nGcells = 6

expn = 10
lown = 10

g = 9.81

dxs = []
L2s = []

expi = 5
ncurr = lown*(2**expi)

sx = -20.0
ex = 20.0
dx = ((ex - sx))/(ncurr-1)
x = arange(sx ,ex+(0.1)*dx,dx)
nx = len(x)

hA,uhA = initialhDB(x,dx,2,1)


dt = 0.5*dx/sqrt(2*g)

hAN,uhAN = Solver(hA,uhA,0,3,g,dx,dt,nGcells,eps)

