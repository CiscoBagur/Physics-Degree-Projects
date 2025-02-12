import numpy as np
import matplotlib.pyplot as plt 

z=[0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,6.0,7.0,8.0,9.0,
10.0,11.0,11.5,12.0,12.5,13.0,13.5,14.0,14.5,15.0,15.5,16.0,16.5,17.0]
B=[0.89,1.20,1.48,1.70,1.89,1.99,2.08,2.13,2.16,2.20,2.22,2.25,2.27,
2.25,2.27,2.25,2.24,2.24,2.23,2.20,2.18,2.21,2.08,2.00,1.90,1.73,1.50,1.23,0.91]
zerror=0.2
Berror=0.02


def B_teorica(x):
	a=2.05
	l=16.4
	N=300
	B_teo=[]
	print(len(x))
	for i in range(len(x)):
		print(x[i])
		valor = 2e-7*(N*100/l)*1000*np.pi*((l-x[i])/np.sqrt(a**2+(l-x[i])**2)+x[i]/np.sqrt(a**2+x[i]**2))
		B_teo.append(valor)
	#B_teo=4e-7*np.pi*((l-x)/sqrt(a**2+(l-x)**2)+x/sqrt(a**2+x**2))
	return B_teo


plt.plot(z,B_teorica(z))
plt.errorbar(z,B,xerr=zerror,yerr=Berror, fmt='.',capsize=3)
plt.ylabel("B(mT)")
plt.xlabel("z(cm)")
plt.title("B vs z Bobina 3 (N=300)")
plt.legend(('B teòrica', 'Mesures experimentals'), loc='lower center')
plt.show()
