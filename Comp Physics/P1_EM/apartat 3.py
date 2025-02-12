import numpy as np
import matplotlib.pyplot as plt 

N = [100,200,300]

quart = [1.72,2.08,2.20]
mitj = [1.90,2.21,2.27]
costat = [1.89,2.18,2.25]

plt.errorbar(N,quart,xerr=0,yerr=0.02, fmt='.',capsize=3,linestyle=":")
plt.errorbar(N,mitj,xerr=0,yerr=0.02, fmt='.',capsize=3,linestyle=":")
plt.errorbar(N,costat,xerr=0,yerr=0.02, fmt='.',capsize=3,linestyle=":")
plt.ylabel("B(mT)")
plt.xlabel("N (voltes)")
plt.title("B vs N")
plt.legend(('$z=\dfrac{1}{4}$', 'centre geomètric', 'al costat del centre'),loc='upper left')
plt.show()