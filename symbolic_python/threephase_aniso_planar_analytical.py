from sympy.solvers import solve
from sympy import Symbol

#In this problem we have
#also imposed a non-zero displacement
#in the y-direction
#Therefore, we should distinguish
#between Ax_alpha which is the total strain in the x-direction
#with Ay_alpha which is the total strain in the y-direction

#Symbol for Length of the system
Lx      = Symbol('Lx')

#Define total strains x-component
Ax_alpha = Symbol('Ax_alpha')
Ax_beta  = Symbol('Ax_beta')
Ax_gamma = Symbol('Ax_gamma')

#Define the intercepts for x-component of displacement
Bx_alpha = Symbol('Bx_alpha')
Bx_beta  = Symbol('Bx_beta')
Bx_gamma = Symbol('Bx_gamma')

#Define the left and right boundary displacements
uxL  = Symbol('uxL')
uxR  = Symbol('uxR')

#Define th interface positions
x1      = Symbol('x1')
x2      = Symbol('x2')

#Define the modulus for alpha phase
#Define (lambda_alpha + 2*mu_alpha)
mod_alpha     = Symbol('mod_alpha')
muprime_alpha = Symbol('muprime_alpha')

#Define the modulus for beta phase
#Define  mod_beta = (lambda_beta + 2*mu_beta)
#Define zeta_beta = 2*lambda_beta + 2*mu_beta + muprime_beta
mod_beta     = Symbol('mod_beta')
muprime_beta = Symbol('muprime_beta')
zeta_beta    = Symbol('zeta_beta')

#Define the modulus for gamma phase
#Define (lambda_gamma + 2*mu_gamma)
mod_gamma     = Symbol('mod_gamma')
muprime_gamma = Symbol('muprime_gamma')

#Symbol for eigenstrains
eT  = Symbol('eT')

#State the system of equations
#Boundary condition at the left end
eq1 = Ax_alpha*(Lx/2) - Bx_alpha

#Continuity of displacement at alpha/beta interface
eq2 = (Ax_alpha*x1 + Bx_alpha) - (Ax_beta*x1 + Bx_beta)

#Continuity of normal stress at alpha/beta interface
eq3 =   mod_alpha * Ax_alpha - mod_beta * Ax_beta \
      + muprime_alpha * Ax_alpha - muprime_beta * Ax_beta \
      + zeta_beta*eT

#Continuity of displacement at beta/gamma interface
eq4 = (Ax_beta*x2 + Bx_beta) - (Ax_gamma*x2 + Bx_gamma)

#Continuity of normal stress at beta/gamma interface
eq5 =   mod_beta * Ax_beta - mod_gamma * Ax_gamma \
      + muprime_beta * Ax_beta - muprime_gamma * Ax_gamma \
      - zeta_beta*eT

eq6 = uxR - Ax_gamma*(Lx/2) - Bx_gamma

sol = solve([eq1, eq2, eq3, eq4, eq5, eq6], [Ax_alpha, Bx_alpha, Ax_beta, Bx_beta, Ax_gamma, Bx_gamma])

#print('Ax_alpha = ', sol[Ax_alpha].factor())
#print('Ax_beta  = ', sol[Ax_beta].factor())
print('Ax_gamma = ', sol[Ax_gamma].factor())

#print('Bx_alpha =', sol[Bx_alpha].factor())
