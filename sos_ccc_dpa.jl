using JuMP
using MosekTools
using DynamicPolynomials
using MultivariatePolynomials
using LinearAlgebra
using TSSOS # important for SOS
# using COSMO

error = 3   # precision digit places
sos_tol = 1 # the maximum degree of unknown SOS polynomials = deg + sos_tol 
a = 1
τ0, τ1,τ2, τ3,τ4, τ5,τ6 = 0.5, a,a, a,a, a,a # S-procedure constants
threshold = 0.01 # for δ

# control closure certificate - DPA
function ccc_dpa(deg,m,finite = true,infinite = false)
    @assert (1 <= m <= N_dpa)
    # synthesize CCC
    # deg: degree of CCC template
    model = Model(optimizer_with_attributes(Mosek.Optimizer)) # COSMO.
    # set_optimizer_attribute(model, MOI.Silent(), true)
    # set_optimizer_attribute(model, "max_iter", 10000) # only for COSMO
    @variable(model, η1 >= threshold) # finite visit
    @variable(model, η2 >= threshold) # infinite visit
    
    # store all functions in list of lists
    allTs, allTcs = [], [] # 3x3
    Ts, Tcs = [], [] # 3
    for i = 1:Q*Q
        T, Tc, Tb = add_poly!(model, vars, deg) # T(x,y), generate polynomial template with given variables and degree
        push!(Ts,T)
        push!(Tcs,Tc)
        if(i%Q == 0)
            push!(allTs, Ts)
            push!(allTcs, Tcs)
            Ts, Tcs = [], []
        end
    end

    # 1st condition for T, >= 0 by default
    for q = 1:Q
        Txy = allTs[q]
        sumTxfi_a, sumTxfi_b = 0, 0
        for i = 1:(m-1)
            Txfi_a = subs(Txy[dpa("a")], y=>f, u=>us_dpa[i]) # T((x,q),(f(x,u_i),q'))
            sumTxfi_a += Txfi_a
            Txfi_b = subs(Txy[dpa("b")], y=>f, u=>us_dpa[i]) # T((x,q),(f(x,u_i),q'))
            sumTxfi_b += Txfi_b
        end
        Txfm_a = subs(Txy[dpa("a")], y=>f, u=>us_dpa[m]) # T((x,q),(f(x,u_m),q'))
        model,_ = add_psatz!(model, Txfm_a + τ0*sumTxfi_a, x, gax_dpa, [], div(maxdegree(Txfm_a + τ0*sumTxfi_a)+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)

        Txfm_b = subs(Txy[dpa("b")], y=>f, u=>us_dpa[m]) # T((x,q),(f(x,u_m),q'))
        model,_ = add_psatz!(model, Txfm_b + τ0*sumTxfi_b, x, gbx_dpa, [], div(maxdegree(Txfm_b + τ0*sumTxfi_b)+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
    end

    # 2nd condition
    for q = 1:Q
        for r = 1:Q
            Txy = allTs[q][r]
            
            Txf_a = subs(allTs[q][dpa("a")],y=>f) # T((x,q),(f(x,u),q'))
            Tfy_a = subs(allTs[dpa("a")][r],x=>f) # T((f(x,u),q'),(y,r))
            model,_ = add_psatz!(model, Txy-τ1*Tfy_a-τ2*Txf_a, [vars;u], [ga_dpa;gin_dpa], [], div(maxdegree(Txy-τ1*Tfy_a-τ2*Txf_a)+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
            
            Txf_b = subs(allTs[q][dpa("b")],y=>f) # T((x,q),(f(x,u),q'))
            Tfy_b = subs(allTs[dpa("b")][r],x=>f) # T((f(x,u),q'),(y,r))
            model,_ = add_psatz!(model, Txy-τ1*Tfy_b-τ2*Txf_b, [vars;u], [gb_dpa;gin_dpa], [], div(maxdegree(Txy-τ1*Tfy_b-τ2*Txf_b)+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
        end
    end

    vzdeg = deg#+1
    # 3rd condition - finite visit
    if (finite)
        Vs, Vcs = [], []
        for q = 1:Q
            Vx, Vxc, Vxb = add_poly!(model, x, vzdeg) # V((x,q)), generate polynomial template V(x,q) with given variables and degree
            if (q in Qvfs) # V(x) ≥ 0
                model,_ = add_psatz!(model, Vx, x, gx_dpa, [], div(deg+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
            end
            push!(Vs, Vx)
            push!(Vcs, Vxc)
        end
        for qvf in Qvfs
            Vy = subs(Vs[qvf],x=>y)
            Vyp = subs(Vs[qvf],x=>yp)
            Tyyp = subs(allTs[qvf][qvf],y=>yp,x=>y)
            for q0 in Q0s
                Txoy = allTs[q0][qvf]
                model,_ = add_psatz!(model, Vy-Vyp-τ3*Txoy-τ4*Tyyp-η1, vars_vf, gp_dpa, [], div(maxdegree(Vy-Vyp-τ3*Txoy-τ4*Tyyp-η1)+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
            end
        end
    end
    
    # 4th condition - infinite visit
    if (infinite)
        Zs, Zcs = [], []
        for q = 1:Q 
            Zx, Zxc, Zxb = add_poly!(model, x, vzdeg) # Z((x,q)), generate polynomial template Z(x,q) with given variables and degree
            if (q in Qninfs) # Z(x) ≥ 0
                model,_ = add_psatz!(model, Zx, x, gx_dpa, [], div(deg+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
            end
            push!(Zs, Zx) # Z(x)'s
            push!(Zcs, Zxc)
        end
        for q0 in Q0s
            for p in Qninfs
                Zy = subs(Zs[p],x=>y) # Z((y,p))
                Txoy = allTs[q0][p] # T((x0,q0),(y,p))
                # Option 1
                for r in Qninfs
                    Zyp = subs(Zs[r],x=>yp) # Z((y',r))
                    Tyyp = subs(allTs[p][r],y=>yp,x=>y) # T((y,p),(y',r))
                    model,_ = add_psatz!(model, Zy-Zyp-τ5*Txoy-τ6*Tyyp-η2, vars_vf, [gninf_dpa;gyp_ninf_dpa], [], div(maxdegree(Zy-Zyp-τ5*Txoy-τ6*Tyyp-η2)+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
                end
                # # Option 2
                # # y ∈ Xvf so L(y) = "b", no need to consider label "a"
                # # Zf_a = subs(Zs[dpa("a")],x=>dyn(y,u)) # Z((f(y,u),p'))
                # # Tyfy_a = subs(allTs[p][dpa("a")],y=>f,x=>y) # T((y,p),(f(y,u),p'))
                # # model,_ = add_psatz!(model, Zy-Zf_a-τ5*Txoy-τ6*Tyfy_a-η2, [vars;u], ?[gninf;gin], [], div(maxdegree(Zy-Zf_a-τ5*Txoy-τ6*Tyfy_a-η2)+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
                # Zf_b = subs(Zs[dpa("b")],x=>dyn(y,u)) # Z((f(y,u),p'))
                # Tyfy_b = subs(allTs[p][dpa("b")],y=>f,x=>y) # T((y,p),(f(y,u),p'))
                # model,_ = add_psatz!(model, Zy-Zf_b-τ5*Txoy-τ6*Tyfy_b-η2, [vars;u], [gninf_dpa;gin_dpa], [], div(maxdegree(Zy-Zf_b-τ5*Txoy-τ6*Tyfy_b-η2)+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
            end
        end
    end


    @objective(model, Min, η1+η2)
    optimize!(model) #solve for coefficients
    status = termination_status(model)  #all_variables(model)
    # monomials for T and V, Z monomials same as V monomials
    Tb = vcat([MultivariatePolynomials.monomials(vars, i) for i = 0:deg]...)
    Vxb = vcat([MultivariatePolynomials.monomials(x, i) for i = 0:vzdeg]...)
    # Zxb = Vxb
    for j = 1:Q # Go through each vector of Tcs
        for k = 1:Q # Go through each element of Tcs
            Tc = value.(allTcs[j][k]) # retrieve value
            for i in eachindex(Tc) # round each constant in Tc
                Tc[i] = round(Tc[i]; digits = error) # round to order of error
            end
            allTcs[j][k] = Tc
            allTs[j][k] = Tc'*Tb # overwrite the polynomials already stored
        end

        if (finite && j in Qvfs)
            Vxc = value.(Vcs[j])  # get the values of each coefficient for V(x)
            for i in eachindex(Vxc)
                Vxc[i] = round(Vxc[i]; digits = error) # round to order of error
            end
            Vcs[j] = Vxc
            Vs[j] = Vxc'*Vxb
        end

        if (infinite && j in Qninfs)
            Zxc = value.(Zcs[j])  # get the values of each coefficient for Z(x)
            for i in eachindex(Zxc)
                Zxc[i] = round(Zxc[i]; digits = error) # round to order of error
            end
            Zcs[j] = Zxc
            Zs[j] = Zxc'*Vxb
        end
    end

    objv = objective_value(model) # primal
    dual = dual_objective_value(model)
    @show value(η1), value(η2), (objv-dual)
    if(abs(objv-dual)<=1e-4 && (string(status) == "SLOW_PROGRESS" || string(status) == "ITERATION_LIMIT"))
        status = "optimal"
    end

    # status might be optimal but if all coefficients approx 10^{-error}, it's essentially 0.
    if (finite && infinite)
        return (objv,status,allTs,Tb,allTcs,Vs,Vxb,Vcs,Zs,Zcs)
    elseif (infinite)
        return (value(η2),status,allTs,Tb,allTcs,Zs,Vxb,Zcs)
    end
    return (value(η1),status,allTs,Tb,allTcs,Vs,Vxb,Vcs) # finite only
end

# Simulation
names = ["limitcycle"]
system_data = ["vars: ","f: ","gx: ","gax: ","gbx: ", "gp: ","gninf: ","gin: "]
max_deg = 3 # does not work past 3 (too big) without enough RAM on device
m = 11 # number of us (inputs) to consider
for name in names
    include("./systems/"*name*".jl"); # load system dynamics info
    file = open("./systems/"*name*"_system_automata.txt", "w"); # open file to write/save system dynamics info
    for (j,k) in zip(system_data, [[vars;u],f,gx_dpa,gax_dpa,gbx_dpa,gp_dpa,gninf_dpa,gin_dpa]) # save these data in readable form
        write(file, j*"{")
        for i = 1:length(k)-1
            write(file, string(k[i])*", ")
        end
        write(file, string(last(k))*"}\n")
    end
    close(file)

    ################ DPA
    println("DPA - PERSISTENCE")
    file = open("./systems/"*name*"_ccc_dpa_persist.txt", "w");
    for deg = 1:max_deg
        stats = @timed data = ccc_dpa(deg, m, true, false) # time execution
        etas, status, Ts,Tmonom,Tcoefs,Vs,Vmonom,Vcoefs = data
        write(file, "poly deg: "*string(deg)*",\tm: "*string(m)*"\n")
        write(file, "status: "*string(status)*",\t η (min): "*string(etas)*"\n")
        for q = 1:Q
            for r = 1:Q
                T, Tc = Ts[q][r], Tcoefs[q][r]
                write(file, "T["*string(q)*string(r)*"]: "*string(T)*"\n")
                write(file, "T_Coef["*string(q)*string(r)*"]: "*string(Tc)*"\n")
            end
        end
        write(file, "Ts Monomials: "*string(Tmonom)*"\n")
        for q in Qvfs
            V, Vc = Vs[q], Vcoefs[q]
            write(file, "V["*string(q)*"]: "*string(V)*"\n")
            write(file, "V_Coef["*string(q)*"]: "*string(Vc)*"\n")
        end
        write(file, "Vs Monomials: "*string(Vmonom)*"\n")
        write(file, "time: "*string(stats.time)*"\n\n")
    end
    close(file)

    ################ DPA
    println("DPA - RECURRENCE")
    file = open("./systems/"*name*"_ccc_dpa_recur.txt", "w");
    for deg = 1:max_deg
        stats = @timed data = ccc_dpa(deg, m, false, true) 
        etas, status, Ts,Tmonom,Tcoefs,Zs,Zmonom,Zcoefs = data
        write(file, "poly deg: "*string(deg)*",\tm: "*string(m)*"\n")
        write(file, "status: "*string(status)*",\t η (min): "*string(etas)*"\n")
        for q = 1:Q
            for r = 1:Q
                T, Tc = Ts[q][r], Tcoefs[q][r]
                write(file, "T["*string(q)*string(r)*"]: "*string(T)*"\n")
                write(file, "T_Coef["*string(q)*string(r)*"]: "*string(Tc)*"\n")
            end
        end
        write(file, "Ts Monomials: "*string(Tmonom)*"\n")
        for q in Qninfs
            Z, Zc = Zs[q], Zcoefs[q]
            write(file, "Z["*string(q)*"]: "*string(Z)*"\n")
            write(file, "Z_Coef["*string(q)*"]: "*string(Zc)*"\n")
        end
        write(file, "Zs Monomials: "*string(Zmonom)*"\n")
        write(file, "time: "*string(stats.time)*"\n\n")
    end
    close(file)

    ################ DPA
    println("DPA - PERSISTENCE + RECURRENCE")
    file = open("./systems/"*name*"_ccc_dpa_persRecur.txt", "w");
    for deg = 1:max_deg
        stats = @timed data = ccc_dpa(deg, m, true, true)
        etas, status, Ts,Tmonom,Tcoefs,Vs,VZmonom,Vcoefs,Zs,Zcoefs = data
        write(file, "poly deg: "*string(deg)*",\tm: "*string(m)*"\n")
        write(file, "status: "*string(status)*",\t Ση (min): "*string(etas)*"\n")
        for q = 1:Q
            for r = 1:Q
                T, Tc = Ts[q][r], Tcoefs[q][r]
                write(file, "T["*string(q)*string(r)*"]: "*string(T)*"\n")
                write(file, "T_Coef["*string(q)*string(r)*"]: "*string(Tc)*"\n")
            end
        end
        write(file, "Ts Monomials: "*string(Tmonom)*"\n")
        for q in Qvfs
            V, Vc = Vs[q], Vcoefs[q]
            write(file, "V["*string(q)*"]: "*string(V)*"\n")
            write(file, "V_Coef["*string(q)*"]: "*string(Vc)*"\n")
        end
        for q in Qninfs
            Z, Zc = Zs[q], Zcoefs[q]
            write(file, "Z["*string(q)*"]: "*string(Z)*"\n")
            write(file, "Z_Coef["*string(q)*"]: "*string(Zc)*"\n")
        end
        write(file, "Vs,Zs Monomials: "*string(VZmonom)*"\n")
        write(file, "time: "*string(stats.time)*"\n\n")
    end
    close(file)
    println("Finished "*name)
end
