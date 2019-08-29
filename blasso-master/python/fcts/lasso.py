import numpy as np


'''
#######################
###   Lasso Fista   ###
#######################
'''


def lasso_primal( y , A , x , lambda_ ):
	return .5*np.linalg.norm(y-A@x)**2+lambda_*np.linalg.norm(x,1)


def lasso_dual_var( y , A , x ):
	return (y-A@x)/np.linalg.norm(A.T.conj()@(y-A@x),np.inf)


def lasso_dual( y , z , lambda_ ):
	return .5*np.linalg.norm(y)**2-.5*np.linalg.norm(y-lambda_*z)**2


def lasso_dual_gap( y , A , x , lambda_):
	z = lasso_dual_var(y , A , x )
	return lasso_primal( y , A , x , lambda_)-lasso_dual( y , z , lambda_ )


def l2_grad( y , A , x ):
	return A.T.conj()@(y-A@x)


def lasso_thsld( x , alpha ):
	return np.maximum((np.abs(x)-alpha),0)*np.exp(1j*np.angle(x))



def fista( y , A , x , lambda_ , L , maxIter=1000 , tol=1.e-6 , disp = True ):
	if(disp):
		print('Fista running...')
	t = 1/(2*L)
	alpha = lambda_*t
	tFISTA = 1
	xFISTA = x
	for iter_ in range(maxIter):
		residual = y-A@xFISTA	
		xold = x
		x = xFISTA+t*A.T.conj()@residual			
		x = lasso_thsld(x,alpha)
		dg = lasso_dual_gap( y , A , x , lambda_)
		if(dg<=tol):
			if(disp):
				print('============')
				print('Fista stopping criteria reached')
				print('Iteration: ',iter_)
				print('Value of the primal function: ',lasso_primal(y,A,x,lambda_))
				print('Value of the dual function: ',lasso_dual(y,lasso_dual_var(y,A,x),lambda_))
				print('Dual gap: ',dg)
				print('============')
			return x
		tFISTA_old = tFISTA
		tFISTA = (1+np.sqrt(1+4*tFISTA**2))/2
		xFISTA = x+(tFISTA_old-1)/tFISTA*(x-xold)
	if(disp):
		print('============')
		print('Maximum iterations reached')
		print('Value of the primal function: ',lasso_primal(y,A,x,lambda_))
		print('Value of the dual function: ',lasso_dual(y,lasso_dual_var(y,A,x),lambda_))
		print('Dual gap: ',dg)
		print('============\n')
	return x

