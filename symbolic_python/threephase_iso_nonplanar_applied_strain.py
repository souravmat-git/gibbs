from sympy.solvers import solve
from sympy import Symbol

#Symbol for outer radius
radius     = Symbol('radius')

#Define total strains
A_alpha = Symbol('A_alpha')
A_beta  = Symbol('A_beta')
A_gamma = Symbol('A_gamma')

#Define the intercepts
B_alpha = Symbol('B_alpha')
B_beta  = Symbol('B_beta')
B_gamma = Symbol('B_gamma')

#Define th interface positions
x1      = Symbol('x1')
x2      = Symbol('x2')

#Define the modulus for alpha phase
#Define (lambda_alpha + 2*mu_alpha)
mod_alpha    = Symbol('mod_alpha')
lambda_alpha = Symbol('lambda_alpha')

#Define the modulus for beta phase
#Define  mod_beta = (lambda_beta + 2*mu_beta)
#Define zeta_beta = 2*lambda_beta + 2*mu_beta + muprime_beta
mod_beta    = Symbol('mod_beta')
lambda_beta = Symbol('lambda_beta')

#Define the modulus for gamma phase
#Define (lambda_gamma + 2*mu_gamma)
mod_gamma    = Symbol('mod_gamma')
lambda_gamma = Symbol('lambda_gamma')

#Define the outer boundary displacements
erg  = Symbol('erg')

#State the system of equations
#Continuity of displacement field at beta/alpha interface
eq1 = (A_alpha - A_beta)*x1 - B_beta/x1

#Continuity of radial stress between alpha and beta phases
eq2 =   (mod_alpha*A_alpha + lambda_alpha*A_alpha) \
      -  mod_beta*(A_beta - B_beta/x1**2) \
      - lambda_beta*(A_beta + B_beta/x1**2)

#Continuity of displacement field at beta/gamma interface
eq3 = (A_beta - A_gamma)*x2 + (B_beta - B_gamma)/x2

#Continuity of radial stress at beta/gamma interface
eq4 =   mod_beta*(A_beta - B_beta/x2**2) \
      + lambda_beta*(A_beta + B_beta/x2**2) \
      - mod_gamma*(A_gamma - B_gamma/x2**2) \
      - lambda_gamma*(A_gamma + B_gamma/x2**2)

#Hoop strain at outer boundary is specified and
#is equal to erg
eq5 = (A_gamma + B_gamma/radius**2) \
      - erg

sol = solve([eq1, eq2, eq3, eq4, eq5], [A_alpha, A_beta, B_beta, A_gamma, B_gamma])

#print('A_alpha = ', sol[A_alpha].factor())
#print('A_beta  = ', sol[A_beta].factor())
print('A_gamma = ', sol[A_gamma].factor())

#print('B_beta =', sol[B_beta].factor())
