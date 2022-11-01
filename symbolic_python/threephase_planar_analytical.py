from sympy.solvers import solve
from sympy import Symbol

#Symbol for Length of the system
Lx      = Symbol('Lx')

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

#Define (2*lambda_alpha + mu_alpha) as mod_alpha
mod_alpha  = Symbol('mod_alpha')

lambda_beta  = Symbol('lambda_beta')
mu_beta      = Symbol('mu_beta')

#Define (2*lambda_beta + mu_beta) as mod_beta
mod_beta      = Symbol('mod_beta')

lambda_gamma = Symbol('lambda_gamma')
mu_gamma     = Symbol('mu_gamma')

#Define (2*lambda_gamma + mu_gamma) as mod_gamma
mod_gamma      = Symbol('mod_gamma')

#Symbol for eigenstrains
eT          = Symbol('eT')

#State the system of equations

#Two equations due to the imposed boundary conditions
eq1 = A_alpha*(Lx/2) - B_alpha
eq6 = - A_gamma*(Lx/2) - B_gamma

#Two equations due to continuity of displacement
eq2 = (A_alpha*x1 + B_alpha) - (A_beta*x1 + B_beta)
eq4 = (A_beta*x2 + B_beta) - (A_gamma*x2 + B_gamma)

#Two equaton due to stress continuity
eq3 = mod_alpha * A_alpha - mod_beta * A_beta \
      + 2 * (lambda_beta + mu_beta)*eT
eq5 = mod_beta * A_beta - 2.0*(lambda_beta + mu_beta)*eT \
       - mod_gamma * A_gamma

sol = solve([eq1, eq2, eq3, eq4, eq5, eq6], [A_alpha, B_alpha, A_beta, B_beta, A_gamma, B_gamma])

#print('A_alpha = ', sol[A_alpha].factor())
#print('A_beta  = ', sol[A_beta].factor())
#print('A_gamma = ', sol[A_gamma].factor())

print('A_alpha =', sol[A_alpha].factor())
