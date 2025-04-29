# NO AUTOMATA
# vector field - system variables
@polyvar x[1:2] y[1:2] yp[1:2] u # yp for persistence
vars = [x;y]
vars_vf = [x;y;yp]

# constants
alpha = 1.5
beta  = 0.5
T = 0.1 # sampling time

# system dynamics - discrete-time Lotka-Volterra type model
function dyn(x,u)
     x1, x2 = x
     f = [x1 + T*(-x2-alpha*x1^2-beta*x1^3),
          x2 + T*(x1-u)]
    return f
end
f = dyn(x,u)

#initial, unsafe sets and state space
# Xo = [0.1,0.6, -0.7,0.7] # initial set (x1,x2)
# Xvf = [-1,0, -1,1] #[-0.9,0.9, -1,-0.8]
# # Xinf = [-3,2.8, -3.3,3.3] # infinite visit
# # c = 0 # epsilon barrier for complement
# # X_Xinf = [Xinf[1]+c,5, -3.3,3.3] # complement of Xinf (X\Xinf)
# X = [-1,1, -1,1] # state space
# U = [-3,3]
#initial, unsafe sets and state space -  jet engine
Xo = [0.1,0.6, -0.7,0.7] # initial set (x1,x2)
Xvf = [-1,0, -1.5,1.5] #[-0.9,0.9, -1,-0.8]
Xinf = [0.5,1, -1.5,1.5] # infinite visit
c = 0 # epsilon barrier for complement
X_Xinf = [-1,Xinf[1]-c, -1.5,1.5] # complement of Xinf (X\Xinf)
X = [-1,1, -1.5,1.5] # state space
U = [-3,3]


# semi-algebraic set description polynomial of the relevant sets, >=0 by default
gx = [x.-X[1:2:length(X)]; X[2:2:length(X)].-x] # state set
g = [gx; 
     y.-X[1:2:length(X)]; X[2:2:length(X)].-y]
# gou = [x.-Xo[1:2:length(Xo)]; Xo[2:2:length(Xo)].-x;  # init (x = xo)
#      y.-Xu[1:2:length(Xu)]; Xu[2:2:length(Xu)].-y] # unsafe (y = xu)
# persistence
gp = [x.-Xo[1:2:length(Xo)]; Xo[2:2:length(Xo)].-x; # init (x = xo)
      y.-Xvf[1:2:length(Xvf)]; Xvf[2:2:length(Xvf)].-y; # finite visit (y, y' = yp)
      yp.-Xvf[1:2:length(Xvf)]; Xvf[2:2:length(Xvf)].-yp]
# Xinf for y
gy_inf = [y.-Xinf[1:2:length(Xinf)]; Xinf[2:2:length(Xinf)].-y]
# complement of Xinf (not inf = ninf)
gy_ninf = [y.-X_Xinf[1:2:length(X_Xinf)]; X_Xinf[2:2:length(X_Xinf)].-y] # not infinite visit (y)
gyp_ninf = [yp.-X_Xinf[1:2:length(X_Xinf)]; X_Xinf[2:2:length(X_Xinf)].-yp] # not infinite visit (y' = yp)
gninf = [x.-Xo[1:2:length(Xo)]; Xo[2:2:length(Xo)].-x; # init (x = xo)
          gy_ninf]
# input bounds
gin = [u.-U[1:2:length(U)]; U[2:2:length(U)].-u] # (old) conservative SOS approach: Pushpak

N = 21
us = LinRange(U[1],U[2],N) # list of inputs. (new)
