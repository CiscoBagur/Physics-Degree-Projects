import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 




x= [5.186e14,5.490e14,6.878e14,7.408e14,8.202e14]
y= [785,901,1504,1689,1921]

yerror = [1,2,1,1,1]
xerror= [0,0,0,0,0]


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

def pesos(yerror):
    triangle= len(y)*sumax2-(sumax)**2
    B=(len(y)*sumaxy-sumax*sumay)/triangle
    w=[]
    for i in range(len(yerror)):
        w_i=1/((yerror[i])**2+(B*xerror[i])**2) #canviar xerror[i] per xerror si la incertesa és fixa
        w.append(w_i)
    return w

def ajustambpesos():
    w= pesos(yerror)
    #sumatori de wi*xi^2
    sumatoriwx2 = 0
    for i in range(len(yerror)):
        sumatoriwx2=sumatoriwx2+w[i]*((x[i])**2)
    #sumatori de w_i * y_i:
    sumatoriwy=0
    for i in range(len(yerror)):
        sumatoriwy= sumatoriwy+(w[i]*y[i])
    #sumatori de w_i * x_i:
    sumatoriwx=0
    for i in range(len(yerror)):
        sumatoriwx=sumatoriwx+(w[i]*x[i])
    #sumatori de w_i * x_i * y_i:
    sumatoriwxy=0
    for i in range(len(yerror)):
        sumatoriwxy=sumatoriwxy+(w[i]*x[i]*y[i])
    #sumatori de w_i
    sumatoriw=0
    for i in range(len(yerror)):
        sumatoriw=sumatoriw+w[i]
    triangle=(sumatoriw*sumatoriwx2)-(sumatoriwx)**2
    A=((sumatoriwx2*sumatoriwy)-(sumatoriwx*sumatoriwxy))/triangle
    sigma_A=np.sqrt(sumatoriwx2/triangle)
    B=((sumatoriw*sumatoriwxy)-(sumatoriwx*sumatoriwy))/triangle
    sigma_B=np.sqrt(sumatoriw/triangle)
    plt.errorbar(x, y, yerr=yerror, xerr=xerror,color="b",fmt='.',capsize=3)
    x_ajust=np.linspace(0,x[len(x)-1])#posar valor mes gran de x per que faja la recta mes llarga
    y_ajust=A+B*x_ajust
    plt.plot(x_ajust,y_ajust,"r",linestyle='dashed')
    plt.ylabel("Vf (mV)")
    plt.xlabel("Freqüència (Hz)")
    plt.title("Vf vs Freqüència")
    plt.ylim(0,2000)
    plt.xlim(0,9e14)
    plt.legend(('Ajust lineal', 'Mesures experimentals'), loc='upper left')
  
    plt.show()
    print("'''''''AJUST AMB INCERTESES VARIABLES (PESOS)'''''''")
    print('A=',A)
    print('B=', B)
    print('sigma_A=', sigma_A)
    print('sigma_B=', sigma_B)
    print("El coeficient de correlació lineal és:", np.corrcoef(x,y)[0][1])

ajustambpesos()