"""
Python code to generate csv files that contain integrals of the basis functions
needed to generate the matrices of the Finite Element Method.

by Jordan Pitt 11/10/2018
"""

"""
############################# IMPORTS #########################################
"""


from sympy import * #Related website: https://www.sympy.org/en/index.html
from IPython.display import display

"""
Function that given a:
    func- function placeholder, just using f
    expandloc - location of desired derivative
    evallocs - location of function evaluations to approximate derivative (relative to expandloc)
    desiredorder - the desired order of the finite difference approximation
    desiredderiv - the desired derivative (first, second, third)
    
    from this we then produce a list
    [a0 a1 a2 ...]
    [f0 f1 f2 ...]
    
    where the finite difference approximation to
    
    d^(desiredderiv) / d(x^desiredderiv) func(x) | x = expandloc  ~ a0*f0 + a1*f1 + a2*f2 + ... 
    
    For this to work we must have:
        desiredderiv + desiredorder = length(evallocs)
"""
def TSFDCoeff(func,expandloc,evallocs,desiredorder,desiredderiv):
    
    if desiredderiv + desiredorder > len(evallocs):
        raise NameError('Not enough evaluation points for FD approximation - require desiredderiv + desiredorder = length(evallocs) ')
    
    elif desiredderiv + desiredorder < len(evallocs):
        raise NameError('Too many evaluation points for FD approximation - require desiredderiv + desiredorder = length(evallocs)')
    
    else:
        #Get TS expansion
        fxaroundexpandloc= func.series(x, x0=expandloc, n=desiredorder +desiredderiv ).doit()
        
        #Get list of derivatives in TS
        derivlist = [func.subs(x,expandloc )] + [diff(func, (x, i)).subs( x, expandloc) for i in range(1,desiredorder +desiredderiv)] 
        
        #Get TS expansions around expandloc for the evallocs
        #fxaroundexpandloc_evals = [ fxaroundexpandloc.subs(x, expandloc + i ).removeO() for i in evallocs ] 
        fxaroundexpandloc_evals = []
        for i in evallocs:
            #For 0, TS is obvious, and also other method doesn't work
            if i  == 0:
                fxaroundexpandloc_evals.append(func.subs(x,expandloc))
            else:
                fxaroundexpandloc_evals.append(fxaroundexpandloc.subs(x, expandloc + i ).removeO())
        
        #Show all locations for function evaluation
        fx_evals = [ func.subs(x, expandloc + i ) for i in evallocs ] 
        
        nds = len(derivlist)
        nexp = len(fxaroundexpandloc_evals )
        i = 0
        
        #Loop over TS expansions for evaluation locations, getting coefficients for
        # derivative terms, by replacing everything by 0, except the desired derived, which is replace with 1
        CoeffMatrix = zeros(nds, nexp)
        for i in range(nexp):
            
            for k in range(nds):
                currcoeff =fxaroundexpandloc_evals[i]
                
                #We replace backwards, since replacing the other way affects derivatives
                for j in range(nds-1,-1,-1):
                    
                    #Replace everything with 0, except when j = k in which case we replace it by 0
                    if j != k:
                        currcoeff = currcoeff.subs(derivlist[j],0)
                    else:
                        currcoeff = currcoeff.subs(derivlist[j],1)
                CoeffMatrix[i,k] = currcoeff
        
        FDCoeffMatrix = CoeffMatrix**(-1)
        Expression = FDCoeffMatrix[desiredderiv,:]*Matrix(fx_evals)
        return Subs(Derivative(func,(x,desiredderiv)),x,expandloc),Expression,FDCoeffMatrix[desiredderiv,:],Matrix(fx_evals)

# symbols
x= symbols('x')
f = Function("f")
h = symbols("dx")
func = f(x)


#Backward
print('Backward Full')
expandloc = symbols('x_{j}')
evallocs = [0,-h, h]
desiredorder = 2
desiredderiv =1
DerivApproximation = TSFDCoeff(func,expandloc,evallocs,desiredorder,desiredderiv)

print('Second order approximation to dq/dx^-_{j+1/2}')
display(DerivApproximation[0])
display(DerivApproximation[1])
print('\n\n\n')

