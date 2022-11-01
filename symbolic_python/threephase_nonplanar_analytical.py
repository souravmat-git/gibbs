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

#Define the left and right boundary displacements
ux_left  = Symbol('ux_left')
ux_right = Symbol('ux_right')

#Define th interface positions
x1      = Symbol('x1')
x2      = Symbol('x2')

#Define the modulus
lambda_alpha = Symbol('lambda_alpha')
mu_alpha     = Symbol('mu_alpha')

#Define
E_alpha      = Symbol('mod_alpha')

lambda_beta  = Symbol('lambda_beta')
mu_beta      = Symbol('mu_beta')

#Define
E_beta      = Symbol('mod_beta')

lambda_gamma = Symbol('lambda_gamma')
mu_gamma     = Symbol('mu_gamma')

#Define
E_gamma      = Symbol('mod_gamma')

#Symbol for eigenstrains
eT          = Symbol('eT')

#State the system of equations

#Continuity of displacement field at beta/alpha interface
eq1 = (A_alpha - A_beta)*x1 - B_beta/x1

#Continuity of radial stress between alpha and beta phases
eq2 =   (E_alpha*A_alpha + lambda_alpha*A_alpha) \
      -  E_beta*(A_beta - B_beta/x1**2) \
      -  lambda_beta*(A_beta + B_beta/x1**2) \
      +  2*(lambda_beta   + mu_beta)*eT

#Continuity of displacement field at beta/gamma interface
eq3 = (A_beta - A_gamma)*x2 + (B_beta - B_gamma)/x2

#Continuity of radial stress at beta/gamma interface
eq4 =   E_beta*(A_beta - B_beta/x2**2) \
      + lambda_beta*(A_beta + B_beta/x2**2) \
      -2*(lambda_beta + mu_beta)*eT\
      - E_gamma*(A_gamma - B_gamma/x2**2) \
      - lambda_gamma*(A_gamma + B_gamma/x2**2)

#Radial stress at outer boundary must be zero
eq5 = E_gamma*(A_gamma - B_gamma/radius**2) \
     + lambda_gamma*(A_gamma + B_gamma/radius**2)

sol = solve([eq1, eq2, eq3, eq4, eq5], [A_alpha, A_beta, B_beta, A_gamma, B_gamma])

#print('A_alpha = ', sol[A_alpha].factor())
#print('A_beta  = ', sol[A_beta].factor())
print('A_gamma = ', sol[A_gamma].factor())

#print('B_beta =', sol[B_beta].factor())
