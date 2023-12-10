import sys
import numpy as np
from scipy.integrate import ode
import scipy.signal
import math
import matplotlib.pyplot as plt
#param
#set up simulator
def f(t, y, arg):
    vin_0 = arg[0]
    ks = arg[1]
    Tin = arg[2]
    c = arg[3]
    A0 = arg[4]
    a = arg[5]
    d = arg[6]
    L = arg[7]
    kk = arg[8]
    Dzero = arg[9]
    n = arg[10]
    X = y[0]
    Y = y[1]
    return [vin_0+A0*math.sin(2*math.pi*t/Tin)-kk*a*X*n*Dzero*(1+a*X/(d+kk))**(n-1)*(1+(40)*Y)**n/((kk+d)*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n)),kk*a*X*n*Dzero*(1+a*X/(d+kk))**(n-1)*(1+40*Y)**n/((kk+d)*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n))-ks*Y]

def jac(t, y, arg):
    vin_0 = arg[0]
    ks = arg[1]
    Tin = arg[2]
    c = arg[3]
    A0 = arg[4]
    a = arg[5]
    d = arg[6]
    L = arg[7]
    kk = arg[8]
    Dzero = arg[9]
    n = arg[10]
    X = y[0]
    Y = y[1]
    return [[-kk*a*n*Dzero*(1+a*X/(d+kk))**(n-1)*(1+(40)*Y)**n/((kk+d)*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n))-kk*a*X*n*Dzero*(n-1)*a*(1+a*X/(d+kk))**(n-2)*(1+(40)*Y)**n/((kk+d)**2*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n))+kk*a*X*n*Dzero*(1+a*X/(d+kk))**(n-1)*(1+(40)*Y)**n*(L*n*(1+a*c*X/d)**(n-1)*c*a/d+n*a*(1+a*X/(d+kk))**(n-1)*(1+40*Y)**n/(kk+d))/((kk+d)*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n)**2),-kk*a*X*n*40*n*Dzero*(1+a*X/(d+kk))**(n-1)*(1+40*Y)**(n-1)/((kk+d)*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n))+kk*a*X*n*Dzero*(1+a*X/(d+kk))**(n-1)*(1+40*Y)**n*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*n*40*(1+40*Y)**(n-1))/((kk+d)*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n)**2)], [kk*a*n*Dzero*(1+a*X/(d+kk))**(n-1)*(1+(40)*Y)**n/((kk+d)*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n))+kk*a*X*n*Dzero*(n-1)*a*(1+a*X/(d+kk))**(n-2)*(1+(40)*Y)**n/((kk+d)**2*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n))-kk*a*X*n*Dzero*(1+a*X/(d+kk))**(n-1)*(1+(40)*Y)**n*(L*n*(1+a*c*X/d)**(n-1)*c*a/d+n*a*(1+a*X/(d+kk))**(n-1)*(1+40*Y)**n/(kk+d))/((kk+d)*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n)**2),kk*a*X*n*40*n*Dzero*(1+a*X/(d+kk))**(n-1)*(1+40*Y)**(n-1)/((kk+d)*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n))-kk*a*X*n*Dzero*(1+a*X/(d+kk))**(n-1)*(1+40*Y)**n*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*n*40*(1+40*Y)**(n-1))/((kk+d)*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n)**2)-ks]]

def main():

	vin_0 = 0.005
	ks = 0.05
	#Tin = 200
#	Tin = int(sys.argv[1])
	c = 0.01
	A0 = float(sys.argv[1])
	a = 2000
	d = 100
	L = 4000000000
	kk = 500
	Dzero = 0.0005
	n = 8

	interval = (10,2000)
	#interval = (10,1000)
	accuracy = 10
	filename = open("/home/purn25/glycolysis_period_A%.4lf.dat"%A0,"w")


	t1 = 100000
	dt = .01
	y0 =[0,0]
	t0 = 0
	TinArray = np.array([])
	YSumArray = np.array([])
	for Tin in np.arange(interval[0],interval[1],accuracy):
	#	params = np.array([vin_0,ks,Tin,c,A0,a,d,L,kk,Dzero,n])
		X = np.zeros( int(t1/dt)+1 )
		Y = np.zeros( int(t1/dt)+1 )
		T = np.zeros( int(t1/dt)+1 )
		r = ode(f, jac).set_integrator('vode', method='bdf', with_jacobian=True)
		r.set_initial_value(y0, t0).set_f_params([vin_0,ks,Tin,c,A0,a,d,L,kk,Dzero,n]).set_jac_params([vin_0,ks,Tin,c,A0,a,d,L,kk,Dzero,n])
		i=0
		while r.successful() and r.t < t1:
			r.integrate(r.t+dt)
			T[i] = r.t
			X[i] = r.y[0]
			Y[i] = r.y[1]
			i=i+1

		Time = []
		for time in range(len(T)):	
			if vin_0+A0*math.sin(2*math.pi*T[time]/Tin)<0:
				Time.append(0)
			else:
				Time.append(vin_0+A0*math.sin(2*math.pi*T[time]/Tin))
		Trange = Time[8000000:]
		Xrange = X[8000000:]
		Yrange = Y[8000000:]
		PeriodList = []
		for peak in range(len(Xrange)):
			if round(Xrange[peak],4) == round(Xrange[0],4)and round(Yrange[peak],4) == round(Yrange[0],4):
#			if Xrange[peak] == Xrange[0] and Yrange[peak] == Yrange[0] and Trange[peak] == Trange[0]:
				PeriodList.append(peak)
#		print(PeriodList)
		Period = []
		perioddiff = np.diff(PeriodList)
		for peakdiff in range(len(perioddiff)):
			if perioddiff[peakdiff] > 2:
				Period.append(perioddiff[peakdiff])
#		print(Period)
		if Period !=[] and len(Period) > 1:
			filename.write("%d\t%lf\t%lf\n"%(Tin,float(np.mean(Period)),float(np.std(Period))))

	filename.close()

	return

if __name__=="__main__":
	main()
