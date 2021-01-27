##### Imports #####
from sympy import *
from sympy.solvers.solveset import linsolve

ga1,ga2,ga3 = symbols('\\gamma_1 \\gamma_2 \\gamma_3')
b1,b2,b3 = symbols('\\beta_1 \\beta_2 \\beta_3')
eps1,eps2,eps3 = symbols('\\epsilon_1 \\epsilon_2 \\epsilon_3')

#We have that ga1 + ga2 + ga3 = 1 and b1,b2,b3>= 0
#We also want either the smallest bi to be largest return of f(bi) and about equal when bi's are equal

#Total value
OT = b1 + b2 + b3

#TTotal Difference
Diff12 = abs(b1 - b2)
Diff13 = abs(b1 - b3)
Diff23 = abs(b2 - b3)
OD = Diff12 + Diff13 + Diff23

OTOD = OT*OD

#Distance From Total (bi - OT)
ODAb1 = b2 + b3
ODAb2 = b1 + b2
ODAb3 = b1 + b3

ODAall = [ODAb1,ODAb2,ODAb3]

#Other Dist From Total |bi - (OT - bi)| = OT - 2*bi
ODANb1 = abs(b1 - (b2 + b3))
ODANb2 = abs(b2 - (b1 + b3))
ODANb3 = abs(b3 - (b2 + b1))

ODANall = [ODANb1,ODANb2,ODANb3]



#Tests:

#all Equal:
OTAllEps1 = OT.subs(b1,eps1).subs(b2,eps1).subs(b3,eps1)
ODAllEps1 = OD.subs(b1,eps1).subs(b2,eps1).subs(b3,eps1)
ODAallEps1 =  [ Exp1.subs(b1,eps1).subs(b2,eps1).subs(b3,eps1) for Exp1 in ODAall] 
OTODE1 = OTOD.subs(b1,eps1).subs(b2,eps1).subs(b3,eps1)
# ODANallEps1 =  [ Exp1.subs(b1,eps1).subs(b2,eps1).subs(b3,eps1) for Exp1 in ODANall] 

#one standout Equal:
OTAllEps12 = OT.subs(b1,eps2).subs(b2,eps1).subs(b3,eps1)
ODAllEps12 = OD.subs(b1,eps2).subs(b2,eps1).subs(b3,eps1)
ODAallEps12 =  [ Exp1.subs(b1,eps2).subs(b2,eps1).subs(b3,eps1) for Exp1 in ODAall] 
OTODE12 = OTOD.subs(b1,eps2).subs(b2,eps1).subs(b3,eps1)
# ODANallEps12 =  [ Exp1.subs(b1,eps2).subs(b2,eps1).subs(b3,eps1) for Exp1 in ODANall] 