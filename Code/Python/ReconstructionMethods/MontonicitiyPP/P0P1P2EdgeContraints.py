from numpy import *
from numpy.linalg import solve,norm
from matplotlib.pyplot import plot,loglog,title,xlabel,ylabel,legend,xlim,ylim



def Monotone(qa):
    n = len(qa)
    IS =  sign(qa[1] - qa[0])

    for i in range(2,n):
        
        CS =sign(qa[i] - qa[i-1]) 
        
        if IS == 0:
            IS = CS
        
        if IS != CS and CS != 0 :
            return False
        
    return True
        
    
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
    if qmin>0 or qmax<0:
        return True
    else:
        return False

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

def PP_IC(x,dx,ldx):
    
    n = len(x)
    q = ones(n)
    
    for i in range(n):
        if (x[i] < 0.1*ldx):
            q[i] = 2 - (x[i])**2
        else:
            q[i] = 1 + 0.01*(x[i])**4
    
    return q  

def EXPPeak_IC(x):
    
    n = len(x)
    q = exp(-abs(x))
    
    return q  

def EXPPeakS_IC(x):
    
    n = len(x)
    # q = exp(-(x**2))
    q = 4 -(x**3)
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
    

def EdgeValuesPjm2j(qajm2,qajm1,qaj):
    return qaj/3 + 5*qajm1/6 - qajm2/6,11*qaj/6 - 7*qajm1/6 + qajm2/3

def EdgeValuesPjm1jp1(qajm1,qaj,qajp1):
    return 5*qaj/6 - qajp1/6 + qajm1/3,5*qaj/6 + qajp1/3 - qajm1/6

def EdgeValuesPjjp2(qaj,qajp1,qajp2):
    return 11*qaj/6 - 7*qajp1/6 + qajp2/3 , qaj/3 + 5*qajp1/6 - qajp2/6

def MaxP2Coeffs(a,b,c):
    return a*(-b/(2*a))**2 +  b*(-b/(2*a))  + c,-b/(2*a)

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
    eps = 10.0**(-12)
    for i in range(nG,n-nG):
        xmh = x[i] - 0.5*dx
        xph = x[i] + 0.5*dx
        
        ndx = (xph - xmh)/np
        xplot = arange(xmh,xph + 0.5*ndx ,ndx)
        xplotm = arange(xmh-dx,xph + 0.5*ndx-dx ,ndx)
        xplotp = arange(xmh+dx,xph + 0.5*ndx+ dx ,ndx)
        
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

        P1jtojp1P =P1jtojp1a11*(xplot - x[i]) + P1jtojp1a10
        P1jm1tojP =P1jm1toja11*(xplot - x[i]) + P1jm1toja10
                    
        Ra1MM = minmodL([P1jtojp1a11,0.5*(P1jtojp1a11 +P1jm1toja11 ),P1jm1toja11])
        Ra0MM = qaj

        RaMMP =  Ra1MM*(xplot - x[i]) + Ra0MM
        
        P0P =  0*(xplot - x[i]) + qaj
        
        
        #Quadratics
        P2jm2toja22,P2jm2toja21,P2jm2toja20,P2jm2tojSI = P2jm2toj(qajm2,qajm1,qaj,dx)   
        P2jm1tojp1a22,P2jm1tojp1a21,P2jm1tojp1a20,P2jm1tojp1SI  = P2jm1tojp1(qajp1,qaj,qajm1,dx)
        P2jtojp2a22,P2jtojp2a21,P2jtojp2a20,P2jtojp2SI =   P2jtojp2(qaj,qajp1,qajp2,dx) 
        
        Pjm2tojLEVj,Pjm2tojREVj =   EdgeValuesPjm2j(qajm2,qajm1,qaj)
        Pjm1tojp1LEVj,Pjm1tojp1REVj = EdgeValuesPjm1jp1(qajm1,qaj,qajp1)
        Pjtojp2LEVj,Pjtojp2REVj =  EdgeValuesPjjp2(qaj,qajp1,qajp2)
        
        MaxLEj = max(Pjm2tojLEVj,Pjm1tojp1LEVj,Pjtojp2LEVj)
        MinLEj = min(Pjm2tojLEVj,Pjm1tojp1LEVj,Pjtojp2LEVj)
        
        MaxREj = max(Pjm2tojREVj,Pjm1tojp1REVj,Pjtojp2REVj)
        MinREj = min(Pjm2tojREVj,Pjm1tojp1REVj,Pjtojp2REVj)

        PbLEVjm1,PbREVjm1 =   EdgeValuesPjm2j(qajm3,qajm2,qajm1)
        PmLEVjm1,PmREVjm1 = EdgeValuesPjm1jp1(qajm2,qajm1,qaj)
        PfLEVjm1,PfREVjm1 =  EdgeValuesPjjp2(qajm1,qaj,qajp1)        

        MaxLEjm1 = max(PbLEVjm1,PmLEVjm1,PfLEVjm1)
        MinLEjm1 = min(PbLEVjm1,PmLEVjm1,PfLEVjm1)
        
        MaxREjm1 = max(PbREVjm1,PmREVjm1,PfREVjm1)
        MinREjm1 = min(PbREVjm1,PmREVjm1,PfREVjm1)

        PbLEVjp1,PbREVjp1 =   EdgeValuesPjm2j(qajm1,qaj,qajp1)
        PmLEVjp1,PmREVjp1 = EdgeValuesPjm1jp1(qaj,qajp1,qajp2)
        PfLEVjp1,PfREVjp1 =  EdgeValuesPjjp2(qajp1,qajp2,qajp3)      

        MaxLEjp1 = max(PbLEVjp1,PmLEVjp1,PfLEVjp1)
        MinLEjp1 = min(PbLEVjp1,PmLEVjp1,PfLEVjp1)
        
        MaxREjp1 = max(PbREVjp1,PmREVjp1,PfREVjp1)
        MinREjp1 = min(PbREVjp1,PmREVjp1,PfREVjp1)

 
    
        Ra2 = 0
        Ra1 = 0
        Ra0 = qaj
        
        print(i,xmh,xph,'Look at Cell')
        print('Left Cell Right Edge Values',PbREVjm1,PmREVjm1,PfREVjm1)
        print('Right Cell LEft Edge Values',PbLEVjp1,PmLEVjp1,PfLEVjp1)
        if (qajm1 >= qaj and qaj >= qajp1):
           
            print('Edge Values j-1/2',MinREjm1,Pjm2tojLEVj,Pjm1tojp1LEVj,Pjtojp2LEVj )
            print('Edge Values j+1/2 ',MaxLEjp1,Pjm2tojREVj,Pjm1tojp1REVj,Pjtojp2REVj  )
            
            if( MinREjm1 >= Pjm2tojLEVj and Pjm2tojREVj >= MaxLEjp1):
                Ra2 = P2jm2toja22
                Ra1 = P2jm2toja21
                Ra0 = P2jm2toja20                
            if( MinREjm1 >= Pjtojp2LEVj and Pjtojp2REVj >= MaxLEjp1):
                Ra2 = P2jtojp2a22
                Ra1 = P2jtojp2a21
                Ra0 = P2jtojp2a20    
            if( MinREjm1 >= Pjm1tojp1LEVj and Pjm1tojp1REVj  >= MaxLEjp1):
                Ra2 = P2jm1tojp1a22
                Ra1 = P2jm1tojp1a21
                Ra0 = P2jm1tojp1a20  
        elif (qajm1 <= qaj and qaj <= qajp1):
            
            print('Edge Values j-1/2',MaxREjm1,Pjm2tojLEVj,Pjm1tojp1LEVj,Pjtojp2LEVj )
            print('Edge Values j+1/2',MinLEjp1,Pjm2tojREVj,Pjm1tojp1REVj,Pjtojp2REVj  )
            if( MaxREjm1 <= Pjm2tojLEVj and Pjm2tojREVj <= MinLEjp1 ):
                Ra2 = P2jm2toja22
                Ra1 = P2jm2toja21
                Ra0 = P2jm2toja20                
            if( MaxREjm1 <= Pjtojp2LEVj and Pjtojp2REVj <= MinLEjp1 ):
                Ra2 = P2jtojp2a22
                Ra1 = P2jtojp2a21
                Ra0 = P2jtojp2a20    
            if( MaxREjm1 <= Pjm1tojp1LEVj and Pjm1tojp1REVj  <= MinLEjp1 ):
                Ra2 = P2jm1tojp1a22
                Ra1 = P2jm1tojp1a21
                Ra0 = P2jm1tojp1a20    
        print('\n\n')
        # #Montone Decreasing - ensure less than average, doesnt work. 
        # qjmh = (qajm1 + qaj)/2
        # qjph = (qaj + qajp1)/2
        # if (qajm1 >= qaj and qaj >= qajp1):
            
        #     if( qjmh >= Pjm2tojLEV and Pjm2tojREV >= qjph):
        #         Ra2 = P2jm2toja22
        #         Ra1 = P2jm2toja21
        #         Ra0 = P2jm2toja20                
        #     if( qjmh >= Pjtojp2LEV and Pjtojp2REV >= qjph):
        #         Ra2 = P2jtojp2a22
        #         Ra1 = P2jtojp2a21
        #         Ra0 = P2jtojp2a20    
        #     if( qjmh >= Pjm1tojp1LEV and Pjm1tojp1REV  >= qjph):
        #         Ra2 = P2jm1tojp1a22
        #         Ra1 = P2jm1tojp1a21
        #         Ra0 = P2jm1tojp1a20  
        # elif (qajm1 <= qaj and qaj <= qajp1):
        #     if( qjmh <= Pjm2tojLEV and Pjm2tojREV <= qjph):
        #         Ra2 = P2jm2toja22
        #         Ra1 = P2jm2toja21
        #         Ra0 = P2jm2toja20                
        #     if( qjmh <= Pjtojp2LEV and Pjtojp2REV <= qjph):
        #         Ra2 = P2jtojp2a22
        #         Ra1 = P2jtojp2a21
        #         Ra0 = P2jtojp2a20    
        #     if( qjmh <= Pjm1tojp1LEV and Pjm1tojp1REV  <= qjph):
        #         Ra2 = P2jm1tojp1a22
        #         Ra1 = P2jm1tojp1a21
        #         Ra0 = P2jm1tojp1a20         

        # print(i,xmh,xph,'Averages', qajm1,qaj,qajp1,  (qajm1 >= qaj and qaj >= qajp1) or (qajm1 <= qaj and qaj <= qajp1) )
        # print(i,xmh,xph,'Edge Values j-2 j',Pjm2tojLEV,Pjm2tojREV, qjmh, qjph, qjmh >= Pjm2tojLEV,Pjm2tojREV >= qjph )
        # print(i,xmh,xph,'Edge Values j-1 j+1',Pjm1tojp1LEV,Pjm1tojp1REV,qjmh, qjph, qjmh >= Pjm1tojp1LEV,Pjm1tojp1REV >= qjph  )
        # print(i,xmh,xph,'Edge Values j j+2',Pjtojp2LEV,Pjtojp2REV, qjmh, qjph, qjmh >= Pjtojp2LEV,Pjtojp2REV >= qjph  )
        
        RaNP2 =  Ra2*(xplot - x[i])**2 + Ra1*(xplot - x[i]) + Ra0
        if i == nG:
            plot(xplot, P0P, '-k',label='Recon P0')  
            
            plot(xplot, P1jm1tojP, '-b',label='Recon P1 -')            
            plot(xplot, P1jtojp1P, '-r',label='Recon P1 +')   
            
            plot(xplot, RaMMP, '-g',label='MinMod P1')   
            
            plot(xplot, RaNP2, '-y',label='Mont P2')   
            
        else:
            plot(xplot, P0P, '-k')  
            
            plot(xplot, P1jm1tojP, '-b')            
            plot(xplot, P1jtojp1P, '-r')   
            
            plot(xplot, RaMMP, '-g')   
            
            plot(xplot, RaNP2, '-y')   





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

# q =  PP_IC(x,dx,dx) 
# qh =  PP_IC(xh,hdx,dx)

plot(xh,qh,'--k', label='Analytic Value')
plot(x,q,'.k', label='Average Values')
PlotAverages(x,q,dx)


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


