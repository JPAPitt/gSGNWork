"""
Python code to calculate analytic expressions of source term in forced solution
for Gaussian with beta values varying in space and time

by Jordan Pitt 11/10/2018
"""

"""
############################# IMPORTS #########################################
"""
from sympy import *
from IPython.display import display

"""
############################# FUNCTIONS #########################################
"""

# Function to write out the analytic solution in a nice format to put in Fortran code
def SympyToCode(Expr,u,h,xdepth):
    NewExpr = Expr
    
    for i in range(xdepth,-1,-1):
        if i == 0 :
            codestru = 'ui'
            codestrh = 'hi'
        elif i == 1 :
            codestru = 'dudx'
            codestrh = 'dhdx'
        else:
            codestru = 'd' + str(i) + 'udx' + str(i)
            codestrh = 'd' + str(i) + 'hdx' + str(i)
        NewExpr= NewExpr.subs(diff(u,(x,i)),codestru ).subs(diff(h,(x,i)),codestrh )
    return NewExpr

# symbols
x,t,ga = symbols('x t ga',real=True)

#Define primitive functions h,u and beta1 and beta2 (since these vary in this forced solution)
h = Function('h')(x,t)
u = Function('u')(x,t)
beta1 = Function('beta1')(x,t)
beta2 = Function('beta2')(x,t)

#Analytic expressions for mass equation (for arbitary h,u,beta1 and beta2)
dhdt = diff(h,t)
fluxh = u*h
dfluxhdx = diff(fluxh,x)
LHSh = dhdt + dfluxhdx  

#Analytic expressions for G equation (transformed conservation of mass equation) (for arbitary h,u,beta1 and beta2)
G = u*h - beta1/2*diff(h**3*diff(u,x),x)

dGdt = diff(G,t)
GReg = -beta2/2*ga*h**2*(h*diff(h,(x,2)) + (diff(h,x)**2)/2) 
dfluxGdx = diff(u*G + ga*h**2/2 - beta1*h**3*diff(u,x)**2 + GReg,x)
LHSG = dGdt + dfluxGdx 

G = u*h - beta1/2*diff(h**3*diff(u,x),x)



#Forced solutions parameters
a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7 = symbols('a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7', positive = True)

#Forced solution terms
phi = x - a2*t
exp1 = exp(-(phi)**2 / (2*a3))

#h expression and derivatives of h
hforce = a0 + a1*exp1
dhfdx = simplify(diff(hforce,x).subs(exp1,'EXPPHI1')).subs(phi,'PHI')
d2hfdx2 = simplify(diff(hforce,(x,2)).subs(exp1,'EXPPHI1')).subs(phi,'PHI')
d3hfdx3 = simplify(diff(hforce,(x,3)).subs(exp1,'EXPPHI1')).subs(phi,'PHI')

#u expression and derivatives of u
uforce = a4*exp1
dufdx = simplify(diff(uforce,x).subs(exp1,'EXPPHI1')).subs(phi,'PHI')
d2ufdx2 = simplify(diff(uforce,(x,2)).subs(exp1,'EXPPHI1')).subs(phi,'PHI')
d3ufdx3 = simplify(diff(uforce,(x,3)).subs(exp1,'EXPPHI1')).subs(phi,'PHI')

#beta expression
beta1force = b1a6*(x - a5*t) + b1a7
beta2force = b2a6*(x - a5*t) + b2a7


#Condensed Form For Coding

#h,u,G expressions and relevant derivatives
Codeh = hforce.subs(exp1,'EXPPHI1').subs(phi,'PHI')
Codeu = uforce.subs(exp1,'EXPPHI1').subs(phi,'PHI')
CodeG = SympyToCode(G,u,h,3)

print('Code For Initial Conditions')
print('h(i) = '+ str(Codeh))
print('u(i) = '+ str(Codeu))
print('dhdx = '+ str(dhfdx))
print('dudx = '+ str(dufdx))
print('d2udx2 = '+ str(d2ufdx2))
print('G(i) = '+ str(CodeG))

print('\n\n\n')
print('Code For Forcing Terms - h')
Codedhdt = simplify((dhdt.subs(h,hforce).doit()).subs(exp1,'EXPPHI1') ).subs(phi,'PHI')
print('dhdt = ' + str(Codedhdt))
print()
Codedfluxhdx = SympyToCode(dfluxhdx,u,h,3)
print('dfluxhdx = ' + str(Codedfluxhdx))
print()

print('Code For Forcing Terms - G')
CodedGdt = dGdt.subs(h,hforce).subs(u,uforce).subs(beta1,beta1force).subs(beta2,beta2force).doit()
CodedGdt = simplify(CodedGdt.subs(uforce,'ui').subs(hforce,'hi') \
            .subs(beta1force,'beta1').subs(beta2force,'beta2')
            .subs(exp1,'EXPPHI1').subs(phi,'PHI')).subs(a2*t - x,'(-PHI)')
CodedGdt = CodedGdt.expand()
print('dGdt = ' + str(CodedGdt))
print()

CodedfluxGdx =dfluxGdx
CodedfluxGdx = SympyToCode(dfluxGdx,u,h,3).subs(beta1,beta1force).subs(beta2,beta2force).doit()
CodedfluxGdx = simplify(CodedfluxGdx.subs(uforce,'ui').subs(hforce,'hi') \
            .subs(beta1force,'beta1').subs(beta2force,'beta2')
            .subs(exp1,'EXPPHI1').subs(phi,'PHI')).subs(a2*t - x,'(-PHI)')
CodedfluxGdx = CodedfluxGdx.expand()
print('dfluxGdx = ' + str(CodedfluxGdx))
print()
print('d2hdx2 = '+ str(d2hfdx2))
print('d3hdx3 = '+ str(d3hfdx3))
print('d3udx3 = '+ str(d3ufdx3))
print('\n\n\n')


