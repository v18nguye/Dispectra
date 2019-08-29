import numpy as np

def Bounds():
	return np.block([[0,1]])

def atoms( theta=.5 , m=51 , sigma=5.e-2 ):
	t = np.linspace(0,1,m)[:,None]
	if(theta.size==1):
		A = np.exp(-(t-theta)**2/(2*sigma**2))
	else:
		theta=theta.flatten('F')[:,None]
		A = np.zeros((m,theta.shape[0]))
		for i in range(theta.shape[0]):
			A[:,i]=np.exp(-(t-theta[i,0])**2/(2*sigma**2))[:,0]
			A[:,i]= A[:,i]/np.linalg.norm(A[:,i])
	return A

def grid( N=100 ):
	return np.linspace( 0 , 1 , N )

def sp( x , y ):
	return np.abs(x.T.conj()@y)

def noise(y,s):
	n = np.random.randn(y.shape[0],y.shape[1])
	return n/np.linalg.norm(n)*np.linalg.norm(y)*10**(-s/20)

def snr(y,n):
	return 20*np.log10(np.linalg.norm(y)/np.linalg.norm(n))	

def signal(k=3,m=51,snr=np.inf):
	theta_y = np.random.rand(1,k)
	x_y = np.random.randn(k,1)
	y = atoms(theta_y,m)@x_y
	n = noise(y,snr)
	y=y+n
	return theta_y,x_y,y
