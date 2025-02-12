import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 

x=[0.00,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00]
y=[0.00,0.23,0.44,0.67,0.91,1.10,1.33,1.54,1.77,1.97,2.21]

yerror = 0.02
xerror= 0.02
#Sumatori de x_i:
sumax=0
for i in range(len(x)):
    sumax=sumax+x[i]


#Sumatori de x_i * y_i:
sumaxy=0
for i in range(len(x)):
    sumaxy=sumaxy+(x[i]*y[i])

#Sumatori de y_i
sumay=0
for i in range(len(x)):
    sumay=sumay+y[i]


#Sumatori de x_i^2
sumax2=0
for i in range(len(x)):
    sumax2=sumax2+(x[i])**2

def ajustsensepesos():
    triangle= len(y)*sumax2-(sumax)**2
    A=(sumax2*sumay-sumax*sumaxy)/triangle
    B=(len(y)*sumaxy-sumax*sumay)/triangle
    erroreq=np.sqrt(yerror**2+(B*xerror)**2)
    sigma_A=erroreq*np.sqrt(sumax2/triangle)
    sigma_B=erroreq*np.sqrt(len(y)/triangle)
    plt.ylabel("$B(mT)$")
    plt.xlabel("$I(A)$")
    plt.title("B vs I Bobina 2 (N=200)")
    plt.ylim(-0.1,2.3)
    plt.xlim(-0.1,1.1)
    x_ajust=np.linspace(0,max(x)+1)
    y_ajust=A+B*x_ajust
    plt.plot(x_ajust,y_ajust,"r",linestyle='dashed')
    plt.errorbar(x,y, yerr=yerror, xerr=xerror,fmt='.',capsize=3)
    plt.legend(('Ajust lineal', 'Mesures experimentals'), loc='upper left')
    plt.show()
    #plt.savefig("grafic.png",dpi=2000)
    print("'''''''AJUST AMB INCERTESES CONSTANTS'''''''")
    print('A=',A)
    print('B=', B)
    print('sigma_A=', sigma_A)
    print('sigma_B=', sigma_B)
    print("El coeficient de correlació lineal és:", (np.corrcoef(x,y)[0][1])**2)

ajustsensepesos()