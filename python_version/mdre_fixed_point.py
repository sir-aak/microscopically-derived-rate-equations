import numpy as np
import mdre_processes as pr

# returns fixed point for given function f, its Jacobian J 
# and an initial guess x using Newton's method
def fixedPointNewton (f, df, x, tol = 1e-10, maxIterations = 100):
	
	xnew = x + 2 * tol
	i    = 0
	
	#~ while np.linalg.norm(np.concatenate((xnew[:1], xnew[2:])) - np.concatenate((x[:1], x[2:]))) > tol and i < maxIterations:
	while np.linalg.norm(xnew - x) > tol and i < maxIterations:
		
		x    = xnew
		#~ xnew = x - np.linalg.inv(df(x, 1e-3)) @ f(x)
		xnew = x - np.linalg.inv(df(x)) @ f(x)
		i   += 1
	
	print("Newton-iterations: ", i)
	
	return (xnew)



# returns fixed point for given function f and an initial guess x 
# using secant method (does not work !!!)
def fixedPointSecant (f, x, tol = 1e-10, maxIterations = 100):
	
	xnew  = x + 2 * tol
	i     = 0
	
	#~ while np.linalg.norm(np.concatenate((xnew[:1], xnew[2:])) - np.concatenate((x[:1], x[2:]))) > tol and i < maxIterations:
	while np.linalg.norm(xnew - x) > tol and i < maxIterations:
		
		xprev = x
		x     = xnew
		xnew  = x - (x - xprev) / (f(x) - f(xprev)) * f(x)	# norm() ?
		i    += 1
	
	print("secant-iterations: ", i)
	
	return (xnew)



# returns fixed point for given function f and an initial guess x 
# using Broyden's method
def fixedPointBroyden (f, df, x, tol = 1e-10, maxIterations = 100):
	
	def updateJacobian (J, f, x, xnew):
		
		dx   = xnew - x
		Jnew = J + np.outer((f(xnew) - f(x) - J @ dx) / np.sum(np.square(dx)), dx)
		
		return (Jnew)
	
	
	J     = df(x)
	xnew  = x + 2 * tol
	i     = 0
	
	#~ while np.linalg.norm(np.concatenate((xnew[:1], xnew[2:])) - np.concatenate((x[:1], x[2:]))) > tol and i < maxIterations:
	while np.linalg.norm(xnew - x) > tol and i < maxIterations:
		
		J    = updateJacobian(J, f, x, xnew)
		x    = xnew
		xnew = x - np.linalg.inv(J) @ f(x)
		i   += 1
	
	print("Broyden-iterations: ", i)
	
	return (xnew)
