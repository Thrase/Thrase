include("read_in.jl")
include("grid.jl")
include("write_out.jl")
include("material_params.jl")
include("operators.jl")
include("solver.jl")

using OrdinaryDiffEq
using Printf
using DelimitedFiles
let

    (T,
     ne,
     Lw,
     r̂,
     l,
     dynamic_flag,
     d_to_s,
     dt_scale,
     ic_file,
     ic_t_file,
     p) = read_params(ARGS[1])

    nn = ne + 1
    Nn = nn^2
    
    # make t_span
    year_seconds = 31556952
    if ic_t_file == "none"
        t_span = (0.0, T * year_seconds)
    else
        f = open(ic_t_file, "r")
        t0 = parse(Float64, readline(f))
        t_span = (t0, T * year_seconds)
    end

    # Get the grid
    μ = 32.038
    ρ = 2.670
    cs = sqrt(μ/ρ)
    xt, yt = transforms_e(Lw, r̂, l)

    @printf "Beginning grid construction...\n"
    metric_time = @elapsed begin
        metrics = create_metrics(ne, μ, xt, yt)
    end
    @printf "Got grid in %f s\n\n" metric_time
    # get fault params
    fc =  metrics.facecoord[2][1]
    
    RS,
    ψ0 = fault_params(fc, μ, cs)
    
    # make output files
    dir_name,
    slip_file,
    slip_rate_file,
    stress_file,
    state_file = make_ss(fc, ne, RS.δNp, ARGS[1])
    station_names = make_stations(dir_name)

    io = (dir_name = dir_name,
          slip_file = slip_file,
          slip_rate_file = slip_rate_file,
          stress_file = stress_file,
          state_file = state_file,
          station_names = station_names,
          pf = [0, 0.0, 0.0])
    
    # Get operators
    R = [-1, 0, 0, 1]
    faces = [0, 2, 3, 4]

    @printf "Beginning operator construction...\n"
    opt_time = @elapsed begin
        ops = operators_dynamic(p, ne, ne, μ, ρ, R, faces, metrics)
    end

    @printf "Got operators in %f seconds\n\n" opt_time
    flush(stdout)

    # initial conditions
    initial_time = @elapsed begin
        # quasi-static
        ψδ = zeros(RS.δNp + nn)
        ψδ[1: RS.δNp] .= ψ0
        vars = (u_prev = zeros(nn^2),
                t_prev = [0.0, 0.0],
                year_seconds = year_seconds,
                Δτ = zeros(nn),
                τ = zeros(nn),
                u = zeros(nn^2),
                v̂_fric = zeros(nn),
                ge = zeros(nn^2))
        
        # dynamic
        q = Array{Float64,1}(undef, 2Nn + 4*nn +  RS.δNp)
    end

    @printf "Got initial conditions in %f seconds\n\n" initial_time
    flush(stdout)
    # set timesteps
    dts = (year_seconds, dt_scale * ops.hmin / (cs))

    # simulation parameters
    params = (reject_step = [false],
              Lw = Lw,
              nn = nn,
              d_to_s = d_to_s,
              vars = vars,
              ops = ops,
              metrics = metrics,
              io = io,
              RS)
              
    # run simulation switching between solvers
    cycle_count = 1

    t_now = 0.0
    @printf "Beginning simulation...\n"
    flush(stdout)
    sim_time = @elapsed begin
        
        #while t_now < t_span[2]
            
            @printf "\tBeginning inter-seismic phase of cycle %d...\n" cycle_count
            flush(stdout)

            inter_time = @elapsed begin

                params.reject_step[1] = false
                stopper = DiscreteCallback(STOPFUN_Q, terminate!)
                prob = ODEProblem(Q_DYNAMIC_SOLVE!, ψδ, t_span, params)
                sol = solve(prob, Tsit5(); isoutofdomain=stepcheck, dt=dts[1],
                            atol = 1e-12, rtol = 1e-12, save_everystep=true,
                            internalnorm=(x, _)->norm(x, Inf), callback=stopper)
                
            end
            
            @printf "\tFinshed inter-seismic phase #%d in %f seconds\n" cycle_count inter_time
            @printf "\tSimulation time: %f years\n" sol.t[end]/year_seconds
        flush(stdout)
            t_now = sol.t[end]
            t_span = (t_now, year_seconds * T)
           
            # update dynamic inital condtions
            q[1 : Nn] .= params.vars.u
        q[Nn + 1 : 2Nn] .= (params.vars.u - params.vars.u_prev)/(sol.t[end] - params.vars.t_prev[1])

        q[2Nn +  1 : 2Nn + nn] = sol.u[end][δNp + 1 : δNp + nn]./2
        
        for 2 in 1:4
            q[2Nn + (i-1)*nn + 1 : 2Nn + i*nn] .= ops.L[i] * params.vars.u
        end
        q[2Nn + 4nn + 1 : 2Nn + 4nn + RS.δNp] .= sol.u[end][1:RS.δNp]

        writedlm("temp_inital.dat", q)
        
            @printf "\tBeginning co-seismic phase of cycle %d...\n" cycle_count
        flush(stdout)
            co_time = @elapsed begin
                t_now = D_timestep!(q, DYNAMIC_SOLVE!, params, dts[2], t_span)
            end
            
            @printf "\tFinished co-seismic phase #%d in %f seconds\n" cycle_count co_time
            @printf "\tSimulation time: %f years\n" t_now/year_seconds
        flush(stdout)
            t_span = (t_now, year_seconds * T)
            vars.t_prev[2] = t_now
            ψδ[1:δNp] .= q[2Nn + 4*nn + 1 : 2Nn + 4*nn + δNp]
            ψδ[δNp + 1 : δNp + nn] .= 2 * q[2Nn +  1 : 2Nn + nn]

            
        #end
    end
    @printf "Ran simulation for %f years and it took %f seconds\n" sol.t[end]/year_seconds sim_time
    flush(stdout)
    nothing
end
