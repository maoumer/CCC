# check barrier from Julia satisfies IBC conditions and if so, plot
# Julia might give wrong answer due to numerical approximation
from z3 import * # SMT solver
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


# system dynamics
def f1(x,u): #jet engine with input
    # constants
    alpha = 1.5
    beta  = 0.5
    T = 0.1; # sampling time
    x1, x2 = x
    return ([x1 + T*(-x2-alpha*x1**2-beta*x1**3),
             x2 + T*(x1-u)])


def f2(x,u): #van der Pol
    # van der Pol
    T = 0.1 # sampling time
    i, v = x
    return ([i + T*v,
            v + T*(-i+u*v*(1-i**2))]) # small 0.1 <= u <= 1.3
    # return ([i + T*(u*(i -1/3*i**3 - v)),
    #          v + T*(1/u*i)]) #  large: 0.5 <= u <= 3


def f3(x,u): # circular limit cycle
    T = 0.1
    x1, x2 = x
    return([x1 + T*(u*x1-x2-x1*(x1**2+x2**2)),
            x2 + T*(x1+u*x2-x2*(x1**2+x2**2))])

### Uncomment/Comment out based on which f is used at the bottom
# #initial, unsafe sets and state space -  jet engine
# Xo = [0.1,0.6, -0.7,0.7] # initial set (x1,x2)
# Xvf = [-1,0, -1.5,1.5] #[-0.9,0.9, -1,-0.8]
# # Xvf = [0,1, -1,0] # unsafe set # keep as Xvf for the rest of code # old: [-2,-1.5, -5,-4]
# Xinf = [0.5,1, -1.5,1.5] # infinite visit
# X = [-1,1, -1.5,1.5] # state space
# U = [-3,3]

# # initial, finite visit sets and state space - van der Pol
# Xo = [3.5,4, 0.5,1] # initial set
# Xvf = [3,3.5, -0.5,0] # finite visit
# Xinf = [-3,2.8, -3.3,3.3] # infinite visit
# X = [-3,4, -3.3,3.3] # state set
# U = [0.1,1.2] # input set

# initial, finite visit sets and state space - circular limit cycle
Xo = [0.8,1, -0.2,0.2] # initial set
Xvf = [0.8,1, 0,0.75] # finite visit
Xinf = [-0.75,0.75, -0.75,0.75] # infinite visit
X = [-0.75,1, -0.75,0.75] # state set
U = [-3,0.5] # input set


x0 = np.array([np.linspace(Xo[0],Xo[1],5), 
                np.linspace(Xo[3],Xo[2],5)])
x02 = np.array([np.linspace(Xo[0],Xo[1],5), 
                np.linspace(Xo[2],Xo[3],5)])
x03 = np.array([np.random.uniform(Xo[0],Xo[1],7),
               np.random.uniform(Xo[2],Xo[3],7)])
x0 = np.concatenate((x0,x02,x03), axis=1)
# print(x0)
x = np.zeros((2,200))
fig, axs = plt.subplots(2,1,sharex=True)
for i in range(x0.shape[1]):
    x[:,0] = x0[:,i]
    for t in range(1,x.shape[1]):
        if (i%3 == 0):
            u = 0.4
        else:
            u = np.random.uniform(U[0],U[1])
        x[:,t] = f3(x[:,t-1],u)

    # t = np.arange(x.shape[1])
    axs[0].stairs(x[0,:],lw=2)
    axs[0].set_ylabel(r"$x_1$")

    axs[1].stairs(x[1,:],lw=2)
    axs[1].set_ylabel(r"$x_2$")
    axs[1].set_xlabel("Time step")

axs[0].grid() 
# axs[0].axis([0,x.shape[1], X[0],X[1]])
axs[0].axhspan(Xo[0], Xo[1], facecolor='#006400', alpha=0.5, label=r"$X_0$")
axs[0].axhspan(Xvf[0], Xvf[1], facecolor='orange',alpha=0.4, label=r"$X_{VF}$") #change to X_VF for lotka-volt
axs[0].axhspan(Xinf[0], Xinf[1], facecolor='gray',alpha=0.5, label=r"$X_{INF}$")

axs[1].grid()
# axs[1].axis([0,x.shape[1], X[2],X[3]])
axs[1].axhspan(Xo[2], Xo[3], facecolor='#006400',alpha=0.5, label=r"$X_0$")  # #90ee90
axs[1].axhspan(Xvf[2], Xvf[3], facecolor='orange',alpha=0.4, label=r"$X_{VF}$") #change to X_VF for lotka-volt
axs[1].axhspan(Xinf[2], Xinf[3], facecolor='gray',alpha=0.5, label=r"$X_{INF}$")
handles, labels = axs[1].get_legend_handles_labels()
fig.legend(handles, labels, loc='center')

savetype = 'svg' # 'png'
# plt.savefig("./media/finitevisit_lv."+savetype, bbox_inches='tight', format=savetype, dpi=1200)
# plt.show()

fig, ax = plt.subplots(1)
for i in range(x0.shape[1]):
    x[:,0] = x0[:,i]
    for t in range(1,x.shape[1]):
        if (i%3 == 0):
            u = 0.4
        else:
            u = np.random.uniform(U[0],U[1])
        x[:,t] = f3(x[:,t-1],u)
    ax.scatter(x[0,:],x[1,:],lw=1)
ax.add_patch(Rectangle((Xo[0], Xo[2]), Xo[1]-Xo[0], Xo[3]-Xo[2],color = '#006400',alpha=0.5,label=r"$X_0$"))
ax.add_patch(Rectangle((Xvf[0], Xvf[2]), Xvf[1]-Xvf[0], Xvf[3]-Xvf[2],color = 'orange',alpha=0.4,label=r"$X_{VF}$")) # change to X_VF for lotka-volt, color to orange
ax.add_patch(Rectangle((Xinf[0], Xinf[2]), Xinf[1]-Xinf[0], Xinf[3]-Xinf[2],color = 'gray',alpha=0.3,label=r"$X_{INF}$")) # change to X_VF for lotka-volt, color to orange
ax.tick_params(axis='both', which='major', labelsize=12)
plt.xlabel(r"$x_1$",fontsize=12)
plt.ylabel(r"$x_2$",fontsize=12)
plt.axis(X)
plt.grid()
plt.legend()
# plt.savefig("./media/trajectories."+savetype, bbox_inches='tight', format=savetype, dpi=1200)
plt.show()
