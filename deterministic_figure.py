import sys
import matplotlib.cm as cm
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
    W = arg[2]
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
    return [vin_0+A0*math.sin(W*t)-kk*a*X*n*Dzero*(1+a*X/(d+kk))**(n-1)*(1+(40)*Y)**n/((kk+d)*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n)),kk*a*X*n*Dzero*(1+a*X/(d+kk))**(n-1)*(1+40*Y)**n/((kk+d)*(L*(1+a*c*X/d)**n+((1+a*X/(d+kk))**n)*(1+40*Y)**n))-ks*Y]

def jac(t, y, arg):
    vin_0 = arg[0]
    ks = arg[1]
    W = arg[2]
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

#	interval = (0.0,10)
	#logscalelist = [0,0.0001,0.000125892,0.0001584893,0.000199526,0.0002511886,0.000316227,0.000398107,0.00050118,0.00063095,0.00079432,0.001,0.00125892,0.001584893,0.00199526,0.002511886,0.00316227,0.00398107,0.0050118,0.0063095,0.0079432,0.01,0.0125892,0.01584893,0.0199526,0.02511886,0.0316227,0.0398107,0.050118,0.063095,0.079432,0.1,0.125892,0.1584893,0.199526,0.2511886,0.316227,0.398107,0.50118,0.63095,0.79432,1]
	logscalelist = [0.00398107]
	#logscalelist = [0.01584893]
	#logscalelist = [0.0050118]
	#logscalelist = [0.0199526]
	#interval = (10,1000)
#	accuracy = 0.05
#fig,biax = plt.subplots()
#fig.set_size_inches(16,9)
	filename = open("/home/purn25/logistic_A%.4lf_W_logY_w0.01995.dat"%A0,"w")


	t1 = 100000
	dt = .01
	y0 =[0,0]
	t0 = 0
	TinArray = np.array([])
	XSumArray = np.array([])
	for W in logscalelist:
	#	print(W)
	#	params = np.array([vin_0,ks,Tin,c,A0,a,d,L,kk,Dzero,n])
		X = np.zeros( int(t1/dt)+1 )
		Y = np.zeros( int(t1/dt)+1 )
		T = np.zeros( int(t1/dt)+1 )
		r = ode(f, jac).set_integrator('vode', method='bdf', with_jacobian=True)
		r.set_initial_value(y0, t0).set_f_params([vin_0,ks,W,c,A0,a,d,L,kk,Dzero,n]).set_jac_params([vin_0,ks,W,c,A0,a,d,L,kk,Dzero,n])
		i=0
		while r.successful() and r.t < t1:
			r.integrate(r.t+dt)
			T[i] = r.t
			X[i] = r.y[0]
			Y[i] = r.y[1]
			i=i+1

	#	Yrange = Y[8000000:]
	#	Ypeaks,_ = scipy.signal.find_peaks(Yrange)
	#	YOpeaks,_ = scipy.signal.find_peaks(-Yrange)
	#	XSum = np.hstack((Yrange[Ypeaks],Yrange[YOpeaks]))
	#	TinArray = np.hstack((TinArray,W*np.ones(len(XSum))))
	#	XSumArray = np.hstack((XSumArray,XSum))
#plot simulation
		Xrange=X[8000000:]
		Yrange=Y[8000000:]
		Xpeaks,_= scipy.signal.find_peaks(Xrange)
		XXpeaks = Xpeaks*dt #normalize time
		XOpeaks,_= scipy.signal.find_peaks(-Xrange)
		Ypeaks,_= scipy.signal.find_peaks(Yrange)
		YOpeaks,_= scipy.signal.find_peaks(-Yrange)
		YYpeaks = Ypeaks*dt #normalize time
#plot Xs and Ys, along with their peaks
##plt.plot(T[99800000:],Y[99800000:]);
		PoincareListX = []
		PoincareListY = []
		XnextStep = []
		YnextStep = []

		for i in range(20):
			PoincareListX.append(X[157746*6+157746*i])
			PoincareListY.append(Y[157746*6+157746*i])
			XnextStep.append(X[157746*7+157746*i])
			YnextStep.append(Y[157746*7+157746*i])
		plt.figure(1)
		plt.subplot(321)
		plt.plot(T[8000000:],Xrange)
		plt.subplot(322)
		plt.plot(T[8000000:],Yrange)
		plt.subplot(323)
		plt.plot(Xrange,Yrange)
		plt.subplot(324)
#plt.scatter(Ypeaks,Yrange[Ypeaks])
#plt.scatter(Xpeaks,Xrange[Xpeaks])
		colors = cm.rainbow(np.linspace(0, 1, 20))
		plt.scatter(PoincareListX, PoincareListY)
		plt.subplot(325)
		plt.scatter(PoincareListX,XnextStep)
		plt.subplot(326)
		plt.scatter(PoincareListY,YnextStep)
		plt.savefig("PoincareXY20_nonchaos_0.0038.pdf")



#	PoincareListX = []
#	PoincareListY = []
#	for j in range(2000):
#		PoincareListX.append(X[100000+314755*i])
#		PoincareListY.append(Y[100000+314755*i])

		
	#for j in range(len(TinArray)):
	#	filename.write("%lf\t%lf\n"%(float(TinArray[j]),float(XSumArray[j])))

	#filename.close()
#	print(Tin*np.ones(len(XSum)),XSum)
#	plt.savefig('/home/purn25/logistic_A%lf'%A0)

	return

if __name__=="__main__":
	main()
