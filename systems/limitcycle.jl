# WITH AND WITHOUT AUTOMATA INFO
# https://math.libretexts.org/Bookshelves/Differential_Equations/A_Second_Course_in_Ordinary_Differential_Equations%3A_Dynamical_Systems_and_Boundary_Value_Problems_(Herman)/03%3A_Nonlinear_Systems/3.08%3A_Limit_Cycles
# vector field - system variables
@polyvar x[1:2] y[1:2] yp[1:2] u # yp for persistence
vars = [x;y]
vars_vf = [x;y;yp] # also for recurrence depending on CCC formulation

T = 0.1; # sampling time
function dyn(x,u) # Example of Hopf bifurcation
    x1, x2 = x
    return([x1 + T*(u*x1-x2-x1*(x1^2+x2^2)),
            x2 + T*(x1+u*x2-x2*(x1^2+x2^2))])
end
f = dyn(x,u)

# no automata
# initial, finite visit sets and state space - circular limit cycle
Xo = [0.8,1, -0.2,0.2] # initial set, [0.8,1]x[-0.2,0.2]
Xvf = [0.8,1, 0,0.75] # finite visit
Xinf = [-0.75,0.75, -0.75,0.75] # infinite visit
c = 0 # epsilon barrier for complement
X_Xinf = [Xinf[2]+c,1, -0.75,0.75] # complement of Xinf (X\Xinf)
X = [-0.75,1, -0.75,0.75] # state set
U = [-3,0.5] # input set


#########################################################
################# NO automata ###########################
#########################################################
# semi-algebraic set description polynomial of the relevant sets, >=0 by default
gx = [x.-X[1:2:length(X)]; X[2:2:length(X)].-x] # state set over x
gy = [y.-X[1:2:length(X)]; X[2:2:length(X)].-y] # state set over y
g = [gx; gy] # state set over x, y

# persistence
gp = [x.-Xo[1:2:length(Xo)]; Xo[2:2:length(Xo)].-x; # init (x = xo)
      y.-Xvf[1:2:length(Xvf)]; Xvf[2:2:length(Xvf)].-y; # finite visit (y, y' = yp)
      yp.-Xvf[1:2:length(Xvf)]; Xvf[2:2:length(Xvf)].-yp]

# recurrence
# Xinf for y
gy_inf = [y.-Xinf[1:2:length(Xinf)]; Xinf[2:2:length(Xinf)].-y]
# complement of Xinf (not inf = ninf)
gy_ninf = [y.-X_Xinf[1:2:length(X_Xinf)]; X_Xinf[2:2:length(X_Xinf)].-y] # not infinite visit (y)
gyp_ninf = [yp.-X_Xinf[1:2:length(X_Xinf)]; X_Xinf[2:2:length(X_Xinf)].-yp] # not infinite visit (y' = yp)
gninf = [x.-Xo[1:2:length(Xo)]; Xo[2:2:length(Xo)].-x; # (x0, y), y not infinite visit
          gy_ninf]

# input bounds
gin = [u.-U[1:2:length(U)]; U[2:2:length(U)].-u]
# discretizing input set
N = 21
us = LinRange(U[1],U[2],N) # list of inputs

#########################################################
################# automata ##############################
#########################################################
# cardinality of parity automata states, initial state, finite visit
Q, Q0s, Qvfs = 3, [1], [2]
Qninfs = [Q0s;Qvfs] # not infinite visit
function dpa(L) # deterministic parity automata transitions
      # for all automata states, same transitions
      if (L == "a") # x in Xinf
            return 3
      else #(L == "b") # x in Xvf
            return 2
      end 
end

# automata
b = 0 # fails if <0 for degree 3 CCC
Xo_dpa = [0.8,1, b,0.2] # initial set
Xvf_dpa = [0.8,1, b,0.75] # finite visit, = not infinite visit
c = 0 # 1e-3
Xinf_dpa = [-0.75,0.8-c, b,0.75] # infinite visit
X_dpa = [-0.75,1, b,0.75] # state set
U = [-3,0.5] # input set

# semialgebraic set for transitions under given labels
gax_dpa = [x.-Xinf_dpa[1:2:length(Xinf_dpa)]; Xinf_dpa[2:2:length(Xinf_dpa)].-x] # letter a, infinite visit gx
gbx_dpa = [x.-Xvf_dpa[1:2:length(Xvf_dpa)]; Xvf_dpa[2:2:length(Xvf_dpa)].-x] # letter b, finite visit gx
gx_dpa = [x.-X_dpa[1:2:length(X_dpa)]; X_dpa[2:2:length(X_dpa)].-x] # all gx over x
gy_dpa = [y.-X_dpa[1:2:length(X_dpa)]; X_dpa[2:2:length(X_dpa)].-y] # all gx over y

# 2nd condition
ga_dpa = [gax_dpa; gy_dpa]
gb_dpa = [gbx_dpa; gy_dpa]

# 3rd condition, persistence
gp_dpa = [x.-Xo_dpa[1:2:length(Xo_dpa)]; Xo_dpa[2:2:length(Xo_dpa)].-x; # init (x = xo)
            y.-Xvf_dpa[1:2:length(Xvf_dpa)]; Xvf_dpa[2:2:length(Xvf_dpa)].-y; # finite visit (y, y' = yp)
            yp.-Xvf_dpa[1:2:length(Xvf_dpa)]; Xvf_dpa[2:2:length(Xvf_dpa)].-yp] 

# 4th condition, gninf = gvf for automata
gninf_dpa = [x.-Xo_dpa[1:2:length(Xo_dpa)]; Xo_dpa[2:2:length(Xo_dpa)].-x;
            y.-Xvf_dpa[1:2:length(Xvf_dpa)]; Xvf_dpa[2:2:length(Xvf_dpa)].-y] #(x0, y)
gyp_ninf_dpa = [yp.-Xvf_dpa[1:2:length(Xvf_dpa)]; Xvf_dpa[2:2:length(Xvf_dpa)].-yp] # not infinite visit (y' = yp)
# input bounds
gin_dpa = [u.-U[1:2:length(U)]; U[2:2:length(U)].-u]

# discretizing input set
N_dpa = 36
us_dpa = LinRange(U[1],U[2],N_dpa) # list of inputs
