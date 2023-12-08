import numpy as np
import sys
import numpy as np
from scipy.integrate import ode
from scipy.integrate import solve_ivp
import scipy.signal
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.spatial.distance import euclidean
#set up simulator
def f(t, y, vin_0,ks,Tin,c,A0,a,d,L,kk,Dzero,n):
	X = y[0]
	Y = y[1]
	if vin_0 + A0*math.sin(2*math.pi*t/Tin)<0:
		inputwave = 0
	else:
		inputwave = vin_0+A0*math.sin(2*math.pi*t/Tin)
	return [inputwave-kk*a*X*n*Dzero*(1+a*X/(d+kk))**(n-1)*(1+(40)*Y)**n/((kk+d)*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n)),kk*a*X*n*Dzero*(1+a*X/(d+kk))**(n-1)*(1+40*Y)**n/((kk+d)*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n))-ks*Y]

def jac(t,y, vin_0,ks,Tin,c,A0,a,d,L,kk,Dzero,n):
    X = y[0]
    Y = y[1]
    return np.array([[-kk*a*n*Dzero*(1+a*X/(d+kk))**(n-1)*(1+(40)*Y)**n/((kk+d)*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n))-kk*a*X*n*Dzero*(n-1)*a*(1+a*X/(d+kk))**(n-2)*(1+(40)*Y)**n/((kk+d)**2*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n))+kk*a*X*n*Dzero*(1+a*X/(d+kk))**(n-1)*(1+(40)*Y)**n*(L*n*(1+a*c*X/d)**(n-1)*c*a/d+n*a*(1+a*X/(d+kk))**(n-1)*(1+40*Y)**n/(kk+d))/((kk+d)*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n)**2),-kk*a*X*n*40*n*Dzero*(1+a*X/(d+kk))**(n-1)*(1+40*Y)**(n-1)/((kk+d)*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n))+kk*a*X*n*Dzero*(1+a*X/(d+kk))**(n-1)*(1+40*Y)**n*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*n*40*(1+40*Y)**(n-1))/((kk+d)*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n)**2)], [kk*a*n*Dzero*(1+a*X/(d+kk))**(n-1)*(1+(40)*Y)**n/((kk+d)*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n))+kk*a*X*n*Dzero*(n-1)*a*(1+a*X/(d+kk))**(n-2)*(1+(40)*Y)**n/((kk+d)**2*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n))-kk*a*X*n*Dzero*(1+a*X/(d+kk))**(n-1)*(1+(40)*Y)**n*(L*n*(1+a*c*X/d)**(n-1)*c*a/d+n*a*(1+a*X/(d+kk))**(n-1)*(1+40*Y)**n/(kk+d))/((kk+d)*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n)**2),kk*a*X*n*40*n*Dzero*(1+a*X/(d+kk))**(n-1)*(1+40*Y)**(n-1)/((kk+d)*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n))-kk*a*X*n*Dzero*(1+a*X/(d+kk))**(n-1)*(1+40*Y)**n*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*n*40*(1+40*Y)**(n-1))/((kk+d)*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n)**2)-ks]])

def func_with_lyaps(t, state, *pars):
    v, U = (state[:2], state[2:].reshape(2,2))
    dv = f(t, v, *pars)
    dU = jac(t, v, *pars) @ U
    return np.hstack((dv,dU.flatten()))

def main():

	vin_0 = 0.005
	ks = 0.05	
	c = 0.01
	A0 = float(sys.argv[1])
	a = 2000
	d = 100
	L = 4000000000
	kk = 500
	Dzero = 0.0005
	n = 8
	interval = (10,2000)
	accuracy = 10
	filename=open('./lyapunov_file/lyapunov_A%.4lf_realsys_Tin.dat'%A0,'w')

# dimension of the system (to keep things general)
	dim_n = 2

# number of Lyapunov exponents
	n_lyap = 2

	assemble = lambda v,U: np.hstack((v,U.flatten()))
	disassemble = lambda state: (state[:dim_n], state[dim_n:].reshape(dim_n,n_lyap))


	for Tin in np.arange(interval[0],interval[1],accuracy):
		v = np.random.random(dim_n)
		U = np.random.random((n_lyap,dim_n))

		lyaps = [] # local Lyapunov exponents

		dt = Tin
		iters = 10000
		tf = iters*dt

		for _ in range(iters):
			sol = solve_ivp(func_with_lyaps,[0, dt],assemble(v,U),method='BDF',t_eval=[dt], args=(vin_0,ks,Tin,c,A0,a,d,L,kk,Dzero,n),max_step =dt)
			v,U = (sol.y.flatten()[:2],sol.y.flatten()[2:].reshape(2,2))
			U,R = np.linalg.qr(U)
			lyaps.append(np.log(abs(R.diagonal()))/dt)

		transient_steps = 200
		result = np.average(lyaps[transient_steps:],axis=0)
		filename.write("%lf\t%lf\t%lf\n"%(Tin,result[0],result[1]))
	filename.close()


	return

if __name__=="__main__":
        main()

