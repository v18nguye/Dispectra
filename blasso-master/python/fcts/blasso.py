import numpy as np
from scipy.optimize import minimize
from . import lasso

'''
#######################
###   Frank-Wolfe   ###
#######################
'''

def blasso_FObj( y , lambda_ , atom , param , coeff , t ):
	return .5*np.linalg.norm(y-atom(param)@coeff)**2+lambda_*t;


def std_FW_update( residual , A , x , t , lambda_ ,t_new=False , Atres_new=False , new_atom=False ):
	if(A.size>1):
		Ax = A@x
	else:
		Ax = np.zeros((residual.shape));
	if(t_new==False):
		c = -Ax;
		gamma = np.min([(np.real(-residual.T.conj()@c)+lambda_*t)/(np.linalg.norm(c)**2),1])
		x = (1-gamma)*x
		t = (1-gamma)*t
	else:
		Adiff_m = -t_new*new_atom*np.exp(1j*np.angle(Atres_new));
        
		c = Adiff_m-Ax
		gamma = np.min([(np.real(-residual.T.conj()@c)+lambda_*(t-t_new))/(np.linalg.norm(c)**2),1])
		if(A.size==1):
			x = -gamma*t_new*np.exp(1j*np.angle(Atres_new))
			t = gamma*t_new
			A = new_atom;
		else:
			x = np.block([[ (1-gamma)*x ],[ -gamma*t_new*np.exp(1j*np.angle(Atres_new)) ]])
			t = (1-gamma)*t+gamma*t_new
			A = np.block([[ A , new_atom ]])
	return A, x, t


def scal_prod( residual , param , atom ):
	return atom(param).T.conj()@residual


def min_scal_prod( residual , param , atom ):
	return -np.sign(np.real(atom(param).T.conj()@residual))*np.real(atom(param).T.conj()@residual)


def FW( y , A_grid , param_grid , atom , B , lambda_ , maxIter=10 , tol=1.e-6 , disp=True):
	if(disp):
		print('Frank-Wolfe running...')
	
	M = np.linalg.norm(y)**2/(2*lambda_);
	t  = 0;
	A = np.array([0])
	param_est=np.array([0])
	x=np.array([0])
	residual = -y

	for iter_ in range(maxIter):
		Atresidual = A_grid.T.conj()@residual
		idx = np.argmax(np.abs(np.real(Atresidual)))
		optim_res = minimize(lambda param_new: min_scal_prod(residual,param_new,atom),(param_grid[idx]),bounds=B[:].T)
		param_new = optim_res.x[0]
		val_new = scal_prod( residual , param_new , atom )

		if(lambda_>=np.abs(val_new)):
			(A , x , t ) = std_FW_update( residual , A , x , t , lambda_ )
		else:
			t_new = M
			new_atom = atom(param_new)
			Atres_new = val_new
			if(A.size==1):
				param_est = param_new[None]
			else:
				param_est = np.block([[ param_est , param_new[None] ]])
			( A , x , t ) = std_FW_update( residual , A , x , t , lambda_ , t_new , Atres_new , new_atom );

		Ax  = atom(param_est)@x
		residual = Ax-y

	if(disp):
		print('============')
		print('Maximum iterations reached')
		print('Value of the blasso primal function: ',blasso_FObj( y , lambda_ , atom , param_est , x , t ))
		print('Value of the lasso primal function: ',blasso_FObj( y , lambda_ , atom , param_est , x , np.sum(np.abs(x)) ))
		print('============\n')
	
	return param_est, x



def SFW( y , A_grid , param_grid , atom , B , lambda_ , maxIter=10 , tol=1.e-6 , disp=True):
	if(disp):
		print('Sliding-Frank-Wolfe running...')
	
	M = np.linalg.norm(y)**2/(2*lambda_);
	t  = 0;
	A = np.array([0])
	param_est=np.array([0])
	x=np.array([0])
	residual = -y

	for iter_ in range(maxIter):
		Atresidual = A_grid.T.conj()@residual
		idx = np.argmax(np.abs(np.real(Atresidual)))
		optim_res = minimize(lambda param_new: min_scal_prod(residual,param_new,atom),(param_grid[idx]),bounds=B[:].T)
		param_new = optim_res.x[0]
		val_new = scal_prod( residual , param_new , atom )

		if(lambda_>=np.abs(val_new)):
			(A , x , t ) = std_FW_update( residual , A , x , t , lambda_ )
		else:
			t_new = M
			new_atom = atom(param_new)
			Atres_new = val_new
			if(A.size==1):
				param_est = param_new[None]
			else:
				param_est = np.block([[ param_est , param_new[None] ]])
			( A , x , t ) = std_FW_update( residual , A , x , t , lambda_ , t_new , Atres_new , new_atom );
		L = np.linalg.eigvals(A.T.conj()@A)
		L = np.max(L)
		x = lasso.fista( y , A , x , lambda_ , L , 100 , tol , False )
		t = np.sum(np.abs(x))
		Ax  = atom(param_est)@x
		residual = Ax-y

	if(disp):
		print('============')
		print('Maximum iterations reached')
		print('Value of the blasso primal function: ',blasso_FObj( y , lambda_ , atom , param_est , x , t ))
		print('Value of the lasso primal function: ',blasso_FObj( y , lambda_ , atom , param_est , x , np.sum(np.abs(x)) ))
		print('============\n')

	return param_est, x
