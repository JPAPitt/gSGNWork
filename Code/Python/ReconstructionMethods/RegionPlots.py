
#Assume we have the expressions at all the necessary points then solved the matrix equation

from numpy import *
from numpy.linalg import solve,norm
from matplotlib.pyplot import plot,loglog



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

def SjplusSjm1(Dj,Djm1,Djp1,alpha):
    Sj = minmodL([alpha*Dj , (Dj + Djp1)/4, alpha*Djp1])
    Sjm1 = minmodL([alpha*Djm1 , (Djm1 + Dj)/4, alpha*Dj])
    return Sj + Sjm1

def TryRanges(RangeDs,alpha):
    n = len(RangeDs)
    for i in range(n):
        for j in range(n):
            SjpSjm1 = SjplusSjm1(RangeDs[i],RangeDs[j],RangeDs[j],alpha)
            
            if (SjpSjm1 <= RangeDs[i]):
                plot(RangeDs[i], RangeDs[j], '.g')  
            else:
                plot(RangeDs[i], RangeDs[j], '.r')  
    

alpha = 1
RangeDs = arange(0,1,0.01)

TryRanges(RangeDs,alpha)