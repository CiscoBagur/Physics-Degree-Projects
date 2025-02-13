import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 

x=[0,0,1,2,3,4,4,5,6,7,7,8,9,10,10,11,12,13,14,15,16,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,32,33,34,35,37,38,39,41,42,43,45,46,48,49,51,52,54,56,58,60,62,63,65,67,70,72,74,77,79,82,85,88,91,94,97,101,105,109,113,117,123,128,135,143,151,159,172,184,198,228,287]
y=[3.801,3.800,3.788,3.776,3.764,3.752,3.739,3.727,3.714,3.701,3.688,3.675,3.662,3.648,3.634,3.620,3.606,3.592,3.577,3.563,3.548,3.532,3.517,3.501,3.485,3.469,3.453,3.436,3.419,3.402,3.384,3.366,3.348,3.330,3.311,3.291,3.272,3.252,3.231,3.210,3.189,3.167,3.145,3.122,3.099,3.075,3.050,3.025,2.999,2.973,2.946,2.918,2.889,2.860,2.829,2.798,2.765,2.732,2.697,2.661,2.624,2.585,2.545,2.503,2.459,2.413,2.365,2.315,2.262,2.206,2.146,2.083,2.016,1.943,1.865,1.781,1.688,1.587,1.473,1.345,1.199,1.027,0.819,0.557,0.200,-0.359,-1.746]

yerror = 0
xerror= 1
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
    plt.xlabel("$t (s)$")
    plt.ylabel("$\ln(\\theta(t)-\\theta_{\infty})$")
    plt.title("$\\theta$ (rad) vs t $(\Delta\Omega = -25\%)$")
    plt.ylim(-2,4)
    plt.xlim(0,300)
    x_ajust=np.linspace(0,max(x)+1)
    y_ajust=A+B*x_ajust
    plt.plot(x_ajust,y_ajust,"r",linestyle='dashed')
    plt.errorbar(x,y, yerr=yerror, xerr=xerror,fmt='.',capsize=3)
    plt.legend(('Ajust lineal', 'Mesures experimentals'), loc='upper right')
    plt.show()
    #plt.savefig("grafic.png",dpi=2000)
    print("'''''''AJUST AMB INCERTESES CONSTANTS'''''''")
    print('A=',A)
    print('B=', B)
    print('sigma_A=', sigma_A)
    print('sigma_B=', sigma_B)
    print("El coeficient de correlació lineal és:", (np.corrcoef(x,y)[0][1])**2)

ajustsensepesos()