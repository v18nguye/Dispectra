import numpy as np

def Bounds():
	return np.block([[-np.pi/2,np.pi/2]])

def atoms( theta=0 , m=51 , Delta=1 , Wavelength=4 ):
	return (1/np.sqrt(m-1))*np.exp(-2*1j*np.pi*Delta/Wavelength*np.arange(-(m-1)/2,(m-1)/2)[:,None]*theta)
	
def grid( N=100 ):
	return np.linspace(-np.pi/2, np.pi/2, N)
		
def sp( x , y ):
	return np.abs(x.T.conj()@y)

def noise(y,s):
	n = np.random.randn(y.shape[0],y.shape[1])
	return n/np.linalg.norm(n)*np.linalg.norm(y)*10**(-s/20)

def snr(y,n):
	return 20*np.log10(np.linalg.norm(y)/np.linalg.norm(n))	

def signal(k=3,m=51,snr=np.inf):
	theta_y = np.pi*np.random.rand(1,k)-np.pi/2
	x_y = np.random.randn(k,1)+1j*np.random.randn(k,1)
	y = atoms(theta_y,m)@x_y
	n = noise(y,snr)
	y=y+n
	return theta_y,x_y,y
	
