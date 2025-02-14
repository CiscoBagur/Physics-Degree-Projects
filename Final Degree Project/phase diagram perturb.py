import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.style.use('seaborn-v0_8-paper')
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
plt.rc('axes', titlesize=15)
plt.rc('axes', labelsize=15)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=12)    # fontsize of the tick labels
plt.rc('ytick', labelsize=12)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize

# Time dependent language shift
def language_shift(t,x):
    return s0*x**a0*(1-x)-(1-s0)*(1-x)**a0*x

def plot_shift(solution,label='Numerical solution'):
    plt.plot(solution.t, solution.y[0], linestyle='--', label=label)
    return

# We've found a perturbative solution arround a=1 for s=1/2
def perturb(t):
	if x0 != 1: return 1/(1+np.exp(-(np.log(x0/(1-x0)))*np.exp((a0-1)*t/2)))
	else: return 1


# Time span parameters
t_span = (0, 1000) 				 	# Simulate for t between 0 and 100
t_eval = np.linspace(*t_span, 5000)  # Time points where we want the solution
x0=0

s0=0.5
ran=100
err_plot = np.eye(ran+1)
a,x = [],[]
for k in range(0,ran+2): x.append(k/ran)
for k in range(0,ran+2): a.append(2*k/ran)

for i in range(0,ran+1):
	for j in range(0,ran+1):
		x0 = i/ran
		a0 = 2*j/ran
		sol = solve_ivp(language_shift, t_span, [x0], t_eval=t_eval)
		error = np.abs(sol.y[0]-perturb(t_eval))
		err_plot[i,j] = np.trapz(error)

err_plot[0,0]=0
err_plot[ran,0]=0
plt.pcolormesh(a,x,err_plot)
plt.xlabel('$\\alpha$')
plt.ylabel('$x_0$')
plt.colorbar()
plt.savefig("diag_bo.png", dpi=500)

