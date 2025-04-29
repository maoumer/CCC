# NO AUTOMATA
# vector field - system variables
@polyvar x[1:2] y[1:2] yp[1:2] u # yp for persistence
vars = [x;y]
vars_vf = [x;y;yp]

T = 0.1; # sampling time

# system dynamics - discrete-time negative-resistance Van der Pol oscillator circuit
function dyn(x,u)
     i, v = x
     return ([i + T*v,
               v + T*(-i+u*v*(1-i^2))]) # small: 0.1 <= u <= 1.3
     # return ([i + T*(u*(i -1/3*i^3 - v)),
     #          v + T*(1/u*i)]) # large: 0.5 <= u <= 3
end
f = dyn(x,u)

# initial, finite visit sets and state space - van der Pol
Xo = [3.5,4, 0.5,1] # initial set
Xvf = [3,3.5, -0.5,0] # finite visit # finite visit
Xinf = [-3,2.8, -3.3,3.3] # infinite visit
c = 0 # epsilon barrier for complement, 0 usually works fine
X_Xinf = [Xinf[2]+c,4, -3.3,3.3] # complement of Xinf (X\Xinf)
X = [-3,4, -3.3,3.3] # state set
U = [0.1,1.2] # input set


# semi-algebraic set description polynomial of the relevant sets, >=0 by default
gx = [x.-X[1:2:length(X)]; X[2:2:length(X)].-x] # state set over x
gy = [y.-X[1:2:length(X)]; X[2:2:length(X)].-y]
g = [gx; gy]

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
gin = [u.-U[1:2:length(U)]; U[2:2:length(U)].-u] # (old) conservative SOS approach

N = 41
us = LinRange(U[1],U[2],N) # list of inputs (new)
