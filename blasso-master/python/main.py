import numpy as np
import matplotlib.pyplot as plt
# import simu.deconv as simu
import simu.doa as simu
import fcts.lasso as lasso
import fcts.blasso as blasso


M = 101
N = 50
K = 5
SNR=np.inf
lambda_lambdaMax = .1

np.random.seed(1)

(theta_y,x_y,y) = simu.signal(K,M,SNR)

theta_A = simu.grid(N)
A = simu.atoms(theta_A,M)

lambdaMax = np.linalg.norm(A.T@y,np.inf)

lambda_ = lambda_lambdaMax*lambdaMax
L = np.linalg.eigvals(A.T.conj()@A)
L = np.max(L)

xinit = np.zeros((N,1))
maxIter = 200
tol = 1.e-6
disp = True


'''
Fista
'''
x_lasso = lasso.fista( y , A , xinit , lambda_ , L , maxIter , tol , disp )


'''
FW-blasso
'''
(param_FW_blasso, x_FW_blasso) = blasso.FW( y , A , theta_A , lambda a:simu.atoms(a,M) , simu.Bounds().T , lambda_ , maxIter , tol , disp )



'''
SFW-blasso
'''
(param_SFW_blasso, x_SFW_blasso) = blasso.SFW( y , A , theta_A , lambda a:simu.atoms(a,M) , simu.Bounds().T , lambda_ , maxIter , tol , disp )


plt.figure()

plt.plot(theta_y.T,np.abs(x_y),'x')
plt.plot(theta_A.T,np.abs(x_lasso))
plt.plot(param_FW_blasso.T,np.abs(x_FW_blasso),'o')
plt.plot(param_SFW_blasso.T,np.abs(x_SFW_blasso),'o')
plt.legend(('GT','fista','FW','SFW'))
plt.show()
