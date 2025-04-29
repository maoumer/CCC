# Computes control closure certificate using SOS for deterministic systems

# include important libraries
using JuMP
using MosekTools
using DynamicPolynomials
using MultivariatePolynomials
using LinearAlgebra
using TSSOS # important for SOS

error = 3   # precision digit places
sos_tol = 1 # the maximum degree of unknown SOS polynomials = deg + sos_tol 
τ0, τ1,τ2, τ3,τ4, τ5,τ6 = 0.5, 1,1, 1,1, 1,1 # S-procedure constants
threshold = 0.1 # for δ

# control closure certificate
function ccc(deg,m,finite = true,infinite = false)
    # deg - degree of CCC template
    # m - # of finite inputs to consider for existential condition
    # finite, infinte - finite or infinite visit conditions
    @assert (1 <= m <= N)
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    # set_optimizer_attribute(model, MOI.Silent(), true)
    @variable(model, η1 >= threshold) # finite visit
    @variable(model, η2 >= threshold) # infinite visit

    # CCC polynomials
    T, Tc, Tb = add_poly!(model, vars, deg) # T(x,y), generate polynomial template T(x,y) with given variables and degree
    
    Txf = subs(T,y=>f) # T(x,f(x,u)) https://juliapackages.com/p/dynamicpolynomials partial substitution
    
    # 1st condition for T, >= 0 by default
    sumTxfi = 0 # -Txf sums fo i = 1 to m-1 using u[i] from finitely many u in U
    for i = 10:10+(m-1)
        Txfi = subs(Txf,u=>us[i]) # T(f(x,u_i),y)
        sumTxfi += Txfi
    end
    Txfm = subs(Txf,u=>us[10+m]) # T(f(x,u_m),y)
    model,_ = add_psatz!(model, Txfm + τ0*sumTxfi, x, gx, [], div(maxdegree(Txfm + τ0*sumTxfi)+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)

    # 2nd condition
    Tfy = subs(T,x=>f) # T(f(x,u),y)
    model,_ = add_psatz!(model, T-τ1*Tfy-τ2*Txf, [vars;u], [g;gin], [], div(maxdegree(T-τ1*Tfy-τ2*Txf)+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
    
    # 3rd condition, T = T(xo,y)
    Tyyp = subs(T,y=>yp,x=>y)
    if (finite) # finite visit
        # persistence
        Vx, Vxc, Vxb = add_poly!(model, x, deg) # generate polynomial template V(x) with given variables and degree
        # V(x) ≥ 0
        model,_ = add_psatz!(model, Vx, x, gx, [], div(deg+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
        Vy = subs(Vx,x=>y) # V(y)
        Vyp = subs(Vx,x=>yp) # V(y')
        model,_ = add_psatz!(model, Vy-Vyp-τ3*T-τ4*Tyyp-η1, vars_vf, gp, [], div(maxdegree(Vy-Vyp-τ3*T-τ4*Tyyp-η1)+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
    end
    # 4th condition, T = T(xo,y)
    if (infinite) # infinite visit
        # recurrence
        Zx, Zxc, Zxb = add_poly!(model, x, deg) # generate polynomial template Z(x) with given variables and degree
        # Z(x) ≥ 0
        model,_ = add_psatz!(model, Zx, x, gx, [], div(deg+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
        Zy = subs(Zx,x=>y) # Z(y)
        # Option 1  - any step
        Zyp = subs(Zx,x=>yp) # Z(y')
        model,_ = add_psatz!(model, Zy-Zyp-τ5*T-τ6*Tyyp-η2, vars_vf, [gninf;gyp_ninf], [], div(maxdegree(Zy-Zyp-τ5*T-τ6*Tyyp-η2)+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
        # # Option 2 - single step version
        # Zyp = subs(Zy,y=>dyn(y,u))
        # Tyyp = subs(Tyyp,yp=>dyn(y,u))
        # model,_ = add_psatz!(model, Zy-Zyp-τ5*T-τ6*Tyyp-η2, [vars;u], [gninf;gin], [], div(maxdegree(Zy-Zyp-τ5*T-τ6*Tyyp-η2)+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
    end

    @objective(model, Min, η1+η2)
    optimize!(model) # solve for coefficients
    status = termination_status(model)
    Tc = value.(Tc)  # get the values of each coefficient for T(x,y)
    for i in eachindex(Tc)
        Tc[i] = round(Tc[i]; digits = error) # round to order of error
    end
    if (finite)
        Vxc = value.(Vxc)  # get the values of each coefficient for V(x)
        for i in eachindex(Vxc)
            Vxc[i] = round(Vxc[i]; digits = error) # round to order of error
        end
    end
    if (infinite)
        Zxc = value.(Zxc)  # get the values of each coefficient for V(x)
        for i in eachindex(Zxc)
            Zxc[i] = round(Zxc[i]; digits = error) # round to order of error
        end
    end

    objv = objective_value(model)
    dual = dual_objective_value(model)
    @show value(η1), value(η2)
    if(abs(objv-dual)<=1e-4 && string(status) == "SLOW_PROGRESS")
        status = "optimal"
    end

    # status might be optimal but if all coefficients approx 10^{-error}, it's essentially 0.
    if (finite && infinite)
        return (objv,status,Tc'*Tb,Tc,Tb,Vxc'*Vxb,Vxc,Vxb,Zxc'*Zxb,Zxc,Zxb)
    elseif (infinite)
        return (value(η2),status,Tc'*Tb,Tc,Tb,Zxc'*Zxb,Zxc,Zxb)
    end
    return (value(η1),status,Tc'*Tb,Tc,Tb,Vxc'*Vxb,Vxc,Vxb) # finite only
end


# Simulation
# names = ["jet_engine","vanderpol"]
names = ["limitcycle"]
system_data = ["vars: ","f: ","g: ","gp: ","gninf: ","gin: "]
max_deg, = 4 # does not work past 4 if not enough RAM on device
m = 6 # number of inputs to consider from finite set
for name in names
    include("./systems/"*name*".jl"); # load system dynamics info
    file = open("./systems/"*name*"_system.txt", "w"); # open file to write/save system dynamics info
    for (j,k) in zip(system_data, [[vars;u],f,g,gp,gninf,gin]) # save these data in readable form
        write(file, j*"{")
        for i = 1:length(k)-1
            write(file, string(k[i])*", ")
        end
        write(file, string(last(k))*"}\n")
    end
    close(file)

    if (name != "vanderpol")
        #### PERSISTENCE - finite visit only
        println("PERSISTENCE") # working with jet engine
        file = open("./systems/"*name*"_ccc_persist.txt", "w");
        for deg = 1:max_deg
            stats = @timed data = ccc(deg, m, true, false) # time execution, (true, false) -> persistence
            eta, status, T, Tcoef, Tmonom, V, Vcoef, Vmonom = data #, Ux
            write(file, "poly deg: "*string(deg)*",\tm: "*string(m)*"\n") #",\tinput deg: "*string(udeg)*
            write(file, "status: "*string(status)*",\t η(min): "*string(eta)*"\n")
            write(file, "T(x,y): "*string(T)*"\n")
            write(file, "Ts Monomials: "*string(Tmonom)*"\n")
            write(file, "Ts Coefficients: "*string(Tcoef)*"\n")
            write(file, "V(x): "*string(V)*"\n")
            write(file, "Vs Monomials: "*string(Vmonom)*"\n")
            write(file, "Vs Coefficients: "*string(Vcoef)*"\n")
            write(file, "time: "*string(stats.time)*"\n\n") 
        end
        close(file)
    end

    if (name != "jet_engine")
        #### RECURRENCE - infinite visit only
        println("RECURRENCE") # working with van der pol
        file = open("./systems/"*name*"_ccc_recur.txt", "w");
        for deg = 1:max_deg-1
            stats = @timed data = ccc(deg, m, false, true) # time execution, (false, true) -> recurrence
            eta, status, T, Tcoef, Tmonom, Z, Zcoef, Zmonom = data
            write(file, "poly deg: "*string(deg)*",\tm: "*string(m)*"\n")
            write(file, "status: "*string(status)*",\t η(min): "*string(eta)*"\n")
            write(file, "T(x,y): "*string(T)*"\n")
            write(file, "Ts Monomials: "*string(Tmonom)*"\n")
            write(file, "Ts Coefficients: "*string(Tcoef)*"\n")
            write(file, "Z(x): "*string(Z)*"\n")
            write(file, "Zs Monomials: "*string(Zmonom)*"\n")
            write(file, "Zs Coefficients: "*string(Zcoef)*"\n")
            write(file, "time: "*string(stats.time)*"\n\n") 
        end
        close(file)
    end

    if (name == "limitcycle")
        #### PERSISTENCE + RECURRENCE - finite and infinite visit
        println("PERSISTENCE + RECURRENCE")
        file = open("./systems/"*name*"_ccc_persRecur.txt", "w");
        for deg = 1:max_deg-1
            stats = @timed data = ccc(deg, m, true, true) # time execution, (true, true) -> persistence + recurrence
            etas, status, T, Tcoef, Tmonom, V, Vcoef, Vmonom, Z, Zcoef, Zmonom = data #, Ux
            write(file, "poly deg: "*string(deg)*",\tm: "*string(m)*"\n") #",\tinput deg: "*string(udeg)*
            write(file, "status: "*string(status)*",\t Ση (min): "*string(etas)*"\n")
            write(file, "T(x,y): "*string(T)*"\n")
            write(file, "Ts Monomials: "*string(Tmonom)*"\n")
            write(file, "Ts Coefficients: "*string(Tcoef)*"\n")
            write(file, "V(x): "*string(V)*"\n")
            write(file, "Vs Monomials: "*string(Vmonom)*"\n")
            write(file, "Vs Coefficients: "*string(Vcoef)*"\n")
            write(file, "Z(x): "*string(Z)*"\n")
            write(file, "Zs Monomials: "*string(Zmonom)*"\n")
            write(file, "Zs Coefficients: "*string(Zcoef)*"\n")
            write(file, "time: "*string(stats.time)*"\n\n") 
        end
        close(file)
    end
    println("Finished "*name)
end
