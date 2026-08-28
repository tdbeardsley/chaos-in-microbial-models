############################################################################################
### Utility Functions ######################################################################
############################################################################################


######################################
### System Specification Functions ###
######################################

function cr_model_co1!(du, u, p, t)
    b = p[1:3]
    k = p[4:6]
    d = p[7:8]
    beta = p[9:12]
    du[1] = (b[1]*u[1]*u[4])/(k[1] + u[4]) - u[1]*d[1]
    du[2] = ((b[2]*u[2]*u[3])/(k[2] + u[3])) + ((b[3]*u[2]*u[4])/(k[3] + u[4])) - d[2]*u[2]
    du[3] = ((beta[1]*u[1]*u[4])/(k[1] + u[4])) - ((beta[2]*u[2]*u[3])/(k[2] + u[3]))
    du[4] = -((beta[3]*u[1]*u[4])/(k[1] + u[4])) - ((beta[4]*u[2]*u[4])/(k[3] + u[4]))
end


function cr_model_co2!(du, u, p, t)
    b1, b2 = p[1:2]
    k1, k2 = p[3:4]
    d1, d2 = p[5:6]
    beta1, beta2, beta3 = p[7:9]
    du[1] = (b1*u[1]*u[4])/(k1 + u[4]) - d1*u[1]
    du[2] = (b2*u[2]*u[3])/(k2 + u[3]) - d2*u[2]
    du[3] = ((beta1*u[1]*u[4])/(k1 + u[4])) - ((beta2*u[2]*u[3])/(k2 + u[3]))
    du[4] = -((beta3*u[1]*u[4])/(k1 + u[4]))
end


function cr_model_co3!(du, u, p, t)
    b1, b2 = p[1:2]
    k1, k2 = p[3:4]
    d1, d2 = p[5:6]
    beta1, beta2, beta3 = p[7:9]
    du[1] = (b1*u[1]*u[4])/(k1 + u[4]) - d1*u[1]
    du[2] = (b2*u[2]*u[3])/(k2 + u[3]) - d2*u[2]
    du[3] = beta1*u[1] - ((beta2*u[2]*u[3])/(k2 + u[3]))
    du[4] = -((beta3*u[1]*u[4])/(k1 + u[4]))
end


function cr_model_mono1!(du, u, p, t)
    b = p[1]
    k = p[2]
    d = p[3]
    beta = p[4]
    du[1] = (b*u[1]*u[2])/(k + u[2]) - u[1]*d
    du[2] = -beta*u[1]*u[2]/(k + u[2])
end


function gLV_sys!(du, u, p, t)
    du[1] = p[1]*u[1] +  p[3]*u[1]*u[1] +  p[4]*u[1]*u[2]
    du[2] = p[2]*u[2] +  p[5]*u[2]*u[1] +  p[6]*u[2]*u[2]
end

###############################################
### Integration functions #####################
###############################################

function dilution_map(x, iv, frac)
    return [frac.*x[1:3]..., frac*x[4] + (1-frac)*iv[4]]
end

function construct_base_forward_dilution_map(prob, p, iv, dt)
    function forward_map(x)
        _prob = DifferentialEquations.remake(prob, u0 = x, p = p, tspan = (0, dt))
        _sol = DifferentialEquations.solve(_prob, Tsit5(), reltol = 1e-16, abstol = 1e-10)
        return dilution_map(_sol.u[end], iv, 0.1)
    end
end

function construct_chaostools_forward_map(system, d_t)

    function chaostools_forward_map!(out, x, p, n)
        prob = DifferentialEquations.ODEProblem(system, x, (0, d_t), p)
        # sol = DifferentialEquations.solve(prob, Feagin14(), saveat = (0, d_t), maxiters = 1e10, reltol = 1e-30, abstol = 1e-20)
        sol = DifferentialEquations.solve(prob, Tsit5(), saveat = (0, d_t), maxiters = 1e8, reltol = 1e-16, abstol = 1e-10)
        out .= dilution_map(sol.u[end], ivs, 0.1)
        return nothing
    end

    return chaostools_forward_map!
end

function integrate_cts_system(
    system,
    parameters,
    data_times,
    t_span,
    initial_values;
    integration_algorithm = Tsit5(),
    reltol = 1e-8,
    abstol = 1e-8,
    maxiters = 1e5
    )

    prob = ODEProblem(
        system,
        initial_values,
        t_span, 
        parameters
    )
    
    sol = solve(
        prob,
        tstops = data_times,
        alg = integration_algorithm,
        reltol = reltol,
        abstol = abstol,
        maxiters = maxiters
    )

    return sol
end


function integrate_dilution_system(
    system,
    parameters,
    initial_values,
    t_span,
    dilution_times,
    dilution_affect;
    integration_algorithm = Tsit5(),
    reltol = 1e-6,
    abstol = 1e-5,
    maxiters = 1e5,
    progress = false,
    save_interval = 20
    )

    cb1 = DifferentialEquations.DiscreteCallback(
        (u, t, integrator) -> t in dilution_times,
        dilution_affect,
        save_positions = (true, true)
    )

    cb2 = PositiveDomain()

    cb = CallbackSet(cb1, cb2)

    prob = DifferentialEquations.ODEProblem(system, initial_values, t_span, parameters)

    save_times = Vector{Float64}
    
    if length(dilution_times) > 1
        save_times = sort!(unique!(vcat(collect(t_span[1]:(dilution_times[2]-dilution_times[1])/save_interval:t_span[end]), dilution_times)))
    else
        save_times = collect(t_span[1]:(t_span[end] - t_span[1])/save_interval:t_span[end])
    end

    sol = solve(
        prob,
        tstops = dilution_times,
        callback = cb,
        alg = integration_algorithm,
        reltol = reltol,
        abstol = abstol,
        maxiters = maxiters,
        progress = progress,
        saveat = save_times
    )

    return sol
end


function dilution_affect_3(fraction, base_ivs)
    function affect!(integrator)
        integrator.u[1:2] = fraction .* integrator.u[1:2]
        integrator.u[3] = fraction * integrator[3]
        integrator.u[4] = (fraction*integrator[4]) + ((1.0-fraction)*base_ivs[4])
    end
    return affect!
end

function get_values_at_times(times, ode_solution_object)
    t_vec = ode_solution_object.t
    val_list = ode_solution_object.u
    time_indices = findall(x -> x in times, t_vec)
    vals = val_list[time_indices]
    return permutedims(reduce(hcat, vals))
end



function get_dilution_values(d_t::Vector, t::Vector, u::Vector; lr = :left)
    d_inds_list = [findall(x -> x in d_t[i], t) for i in eachindex(d_t)]
    if lr == :left
        left_inds = [d_inds_list[i][1] for i in eachindex(d_inds_list)]
        return u[left_inds]
    else
        right_inds = [d_inds_list[i][end] for i in eachindex(d_inds_list)]
        return u[right_inds]
    end
end


function get_dilution_values(d_t::Vector, t::Vector, u::Matrix; lr = :left)
    d_inds = findall(x -> x in d_t, t)
    if lr == :left
        return u[d_inds, :][1:2:end-1, :]
    else
        return u[d_inds, :][2:2:end, :]
    end
end


function get_dilution_values(d_t::Vector, s::ODESolution; lr = :left)
    t = s.t
    u = s.u
    return get_dilution_values(d_t, t, u, lr = lr)
end


function cl(x::Vector)
    return x ./ sum(x)
end


function glv_forward_map(problem)
    function forward_map(x, dt, p)
        prob = remake(problem, u0 = [x, 1.0 - x], tspan = (0, dt), p = p)
        sol = DifferentialEquations.solve(prob , Tsit5(), reltol = 1e-16, abstol = 1e-10)
        return cl(sol.u[end])[1]
    end
    return forward_map
end

function glv_chaostools_forward_map(problem, dt)
    function forward_map!(out, x, p, t)
        prob = remake(problem, u0 = [x[1], 1.0-x[1]], tspan = (0, dt), p = p)
        sol = DifferentialEquations.solve(prob, Tsit5(), reltol = 1e-16, abstol = 1e-10)
        out .= [(cl(sol.u[end]))[1]]
        return nothing
    end
    return forward_map!
end


### Example Loss Function 
function loss_maker_cts_system_fit_co(p, system, data, data_times, tspan, ivs, w)
    
    vals = integrate_cts_system(
        system,
        p, 
        data_times,
        tspan,
        ivs,
        integration_algorithm = Tsit5(),
        reltol = 1e-6,
        abstol = 1e-5,
        maxiters = 1e7
    )

    pred_vals = get_values_at_times(data_times, vals)

    if Int(vals.retcode) == 1
        diffs = (sqrt.(w)) .* (log.(data) .- log.(pred_vals[:, 1:2]))
        # diffs = (sqrt.(w)) .* (data .- pred_vals[:, 1:2])
        return sqrt(sum(abs2, diffs))
    else
        return Inf
    end

end



















