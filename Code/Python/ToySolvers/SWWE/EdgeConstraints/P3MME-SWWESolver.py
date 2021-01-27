
#Assume we have the expressions at all the necessary points then solved the matrix equation

from numpy import *
from numpy.linalg import solve,norm
from matplotlib.pyplot import plot,loglog

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
    

def minmodL(q):
    qmin = min(q)
    qmax = max(q)
    
    # min_index = q.index(qmin)
    # max_index = q.index(qmax)
    if qmin>0:
        return qmin
    if qmax<0 :
        return qmax
    else:
        return 0

def P1jtojp1(qaj,qajp1,dx):
    return (qajp1 - qaj)/ dx ,qaj,(qaj**2 - 2*qaj*qajp1 + qajp1**2)

def P1jm1toj(qaj,qajm1,dx):
    return (qaj - qajm1)/ dx ,qaj,(qaj**2 - 2*qaj*qajm1 + qajm1**2)


def EdgeValuesPjm2j(qajm2,qajm1,qaj):
    return qaj/3 + 5*qajm1/6 - qajm2/6,11*qaj/6 - 7*qajm1/6 + qajm2/3

def EdgeValuesPjm1jp1(qajm1,qaj,qajp1):
    return 5*qaj/6 - qajp1/6 + qajm1/3,5*qaj/6 + qajp1/3 - qajm1/6

def EdgeValuesPjjp2(qaj,qajp1,qajp2):
    return 11*qaj/6 - 7*qajp1/6 + qajp2/3 , qaj/3 + 5*qajp1/6 - qajp2/6

def P2AvgEdges(qaj,qjmh,qjph,dx):
    pa =3*(qjmh + qjph - 2*qaj)/ (dx**2)
    pb = (qjph - qjmh)/dx
    pc = (6*qaj - qjmh - qjph)/(4)
    return pa,pb,pc

def EvalPoly(Coeffs,x):
    PolyVal = Coeffs[0]
    n = len(Coeffs)
    for i in range(1,n):
        PolyVal = PolyVal + Coeffs[i]*x**i
    return PolyVal

def ReconjBE(qajm3,qajm2,qajm1,qaj,qajp1,qajp2,qajp3,dx,eps,j):
    
    # P2jm2toja22,P2jm2toja21,P2jm2toja20,P2jm2tojSI = P2jm2toj(qajm2,qajm1,qaj,dx)   
    # P2jm1tojp1a22,P2jm1tojp1a21,P2jm1tojp1a20,P2jm1tojp1SI  = P2jm1tojp1(qajp1,qaj,qajm1,dx)
    # P2jtojp2a22,P2jtojp2a21,P2jtojp2a20,P2jtojp2SI =   P2jtojp2(qaj,qajp1,qajp2,dx) 
    
    Pjm2tojLEVj,Pjm2tojREVj =   EdgeValuesPjm2j(qajm2,qajm1,qaj)
    Pjm1tojp1LEVj,Pjm1tojp1REVj = EdgeValuesPjm1jp1(qajm1,qaj,qajp1)
    Pjtojp2LEVj,Pjtojp2REVj =  EdgeValuesPjjp2(qaj,qajp1,qajp2)
    
    P1jtojp1a11,P1jtojp1a10,P1jtojp1SI  =  P1jtojp1(qaj,qajp1,dx)
    P1jm1toja11,P1jm1toja10,P1jm1tojSI  =  P1jm1toj(qaj,qajm1,dx)
                
    Ra1MM = minmodL([P1jtojp1a11,0.5*(P1jtojp1a11 +P1jm1toja11 ),P1jm1toja11])
    
    qmh = Ra1MM*(-dx/2) + qaj
    qph = Ra1MM*(dx/2) + qaj
    
    RELV = minmodL([Pjm2tojLEVj - qmh,Pjm1tojp1LEVj - qmh,Pjtojp2LEVj - qmh]) + qmh
    RERV = minmodL([Pjm2tojREVj - qph,Pjm1tojp1REVj - qph,Pjtojp2REVj - qph]) + qph
 
    
    #Cubics
    P3jm3toja33,P3jm3toja32,P3jm3toja31,P3jm3toja30,P3jm3tojSI = P3jm3toj(qaj,qajm1,qajm2,qajm3,dx)
    P3jm2tojp1a33,P3jm2tojp1a32,P3jm2tojp1a31,P3jm2tojp1a30,P3jm2tojp1SI = P3jm1tojp2(qaj,qajm1,qajp1,qajp2,dx)
    P3jm1tojp2a33,P3jm1tojp2a32,P3jm1tojp2a31,P3jm1tojp2a30,P3jm1tojp2SI = P3jm2tojp1(qaj,qajm1,qajm2,qajp1,dx)
    P3jtojp3a33,P3jtojp3a32,P3jtojp3a31,P3jtojp3a30,P3jtojp3SI = P3jtojp3(qaj,qajp1,qajp2,qajp3,dx)
    
    
    Pjm3tojLEV = EvalPoly([P3jm3toja30,P3jm3toja31,P3jm3toja32,P3jm3toja33],-dx/2)
    Pjm3tojREV = EvalPoly([P3jm3toja30,P3jm3toja31,P3jm3toja32,P3jm3toja33],dx/2)
    
    Pjm2tojp1LEV = EvalPoly([P3jm2tojp1a30,P3jm2tojp1a31,P3jm2tojp1a32,P3jm2tojp1a33],-dx/2)
    Pjm2tojp1REV = EvalPoly([P3jm2tojp1a30,P3jm2tojp1a31,P3jm2tojp1a32,P3jm2tojp1a33],dx/2)        
    
    Pjm1tojp2LEV = EvalPoly([P3jm1tojp2a30,P3jm1tojp2a31,P3jm1tojp2a32,P3jm1tojp2a33],-dx/2)
    Pjm1tojp2REV = EvalPoly([P3jm1tojp2a30,P3jm1tojp2a31,P3jm1tojp2a32,P3jm1tojp2a33],dx/2)     
    
    Pjtojp3LEV = EvalPoly([P3jtojp3a30,P3jtojp3a31,P3jtojp3a32,P3jtojp3a33],-dx/2)
    Pjtojp3REV = EvalPoly([P3jtojp3a30,P3jtojp3a31,P3jtojp3a32,P3jtojp3a33],dx/2)   
     
    REL3V = minmodL([Pjm3tojLEV - RELV,Pjm1tojp2LEV - RELV,Pjm2tojp1LEV - RELV,Pjtojp3LEV - RELV]) + RELV
    RER3V = minmodL([Pjm3tojREV - RERV,Pjm1tojp2REV - RERV,Pjm2tojp1REV - RERV,Pjtojp3REV - RERV]) + RERV
    return REL3V,RER3V

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
    
    if hir < hip1l:
        print(j,hir,hi1pl,ha[j-1],haj[j],ha[j+1])
    
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

        # if hir < hip1l:
        #     print(j,hir,hip1l,ha[j-1],ha[j],ha[j+1])  
        #     print(foh,fih,fouh,fiuh)
            
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
    
    hAppp = hApp/4 + 3*hA/4
    uhAppp = uhApp/4 + 3*uhA/4
    
    hApppp,uhApppp = EulerStep(hAppp,uhAppp,g,dx,nGcells,dt,eps)
    
    hAn = hA/3 + 2*hApppp/3
    uhAn = uhA/3 + 2*uhApppp/3
    
    return hAn ,uhAn
    

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

