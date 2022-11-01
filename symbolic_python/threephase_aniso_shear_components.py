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

#Define total strains y-component
Ay_alpha = Symbol('Ay_alpha')
Ay_beta  = Symbol('Ay_beta')
Ay_gamma = Symbol('Ay_gamma')

#Define the intercepts for y-component of displacement
By_alpha = Symbol('By_alpha')
By_beta  = Symbol('By_beta')
By_gamma = Symbol('By_gamma')

#Define the left and right boundary displacements
uyL  = Symbol('uyL')
uyR  = Symbol('uyR')

#Define th interface positions
x1      = Symbol('x1')
x2      = Symbol('x2')

#Define the modulus for alpha phase
#Define (lambda_alpha + 2*mu_alpha)
mu_alpha  = Symbol('mu_alpha')
mu_beta   = Symbol('mu_beta')
mu_gamma  = Symbol('mu_gamma')


#State the system of equations
#Boundary condition at the left end
eq1 = Ay_alpha*(Lx/2) - By_alpha

#Continuity of displacement at alpha/beta interface
eq2 = (Ay_alpha*x1 + By_alpha) - (Ay_beta*x1 + By_beta)

#Continuity of normal stress at alpha/beta interface
eq3 =  2 * mu_alpha * Ay_alpha -  2 * mu_beta * Ay_beta

#Continuity of displacement at beta/gamma interface
eq4 = (Ay_beta*x2 + By_beta) - (Ay_gamma*x2 + By_gamma)

#Continuity of normal stress at beta/gamma interface
eq5 = 2 * mu_beta * Ay_beta - 2 * mu_gamma * Ay_gamma

eq6 = uyR - Ay_gamma*(Lx/2) - By_gamma

sol = solve([eq1, eq2, eq3, eq4, eq5, eq6], [Ay_alpha, By_alpha, Ay_beta, By_beta, Ay_gamma, By_gamma])

#print('Ay_alpha = ', sol[Ay_alpha].factor())
print('Ay_beta  = ', sol[Ay_beta].factor())
#print('Ay_gamma = ', sol[Ay_gamma].factor())

#print('Bx_alpha =', sol[Bx_alpha].factor())
