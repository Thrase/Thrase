include("odefun.jl")
using Plots
using SparseArrays
using LinearAlgebra
using DelimitedFiles
using Interpolations

include("ops_stripped.jl")
 

function main()

    # number of grid points in each dimension
    Ny = 160
    Nz = 160

    stride_space = 1 # write-out every stride_space grid points along fault
    stride_time = 5  # write-out every stride_time time steps

    # Physical domain: (y, z) in (0, Ly) x (0, Lz)
    Ly = 80  # (km)
    Lz = 80  # (km)

    sim_years = 300.  # (years) to simulate.

    # rate-and-state and elasticity parameters:
    Vp = 1e-9   # (m/s)
    ρ = 2.670   # (1000 kg/m^3 - check this)
    cs = 3.464  # (m/s)
    μ = cs^2 * ρ  # (GPa)
    σn = 50       # (MPa)
    RSamin = 0.01
    RSamax = 0.025
    RSb = 0.015
    RSDc = 0.032  # (m)
    RSf0 = 0.6
    RSV0 = 1e-6    # (m/s)
    RSVinit = 1e-9 # (m/s)
    RSH1 = 15      # (km)
    RSH2 = 18      # (km)
    RSWf = 40      # (km)
    μshear = cs^2 * ρ  # (GPa)
    η = μshear / (2 * cs)  # TODO: fill me in



    # SBP finite difference interior order of accuracy:
    SBPp   = 2


    # Physical Domain: (y, z) in (0, Ly) x (0, Lz)
    y = Array(LinRange(0, Ly, Ny+1))  
    z = Array(LinRange(0, Lz, Nz+1))

    # create finite difference operators on computational domain:
    (M̃, F, τ, H̃, HfI_FT) = get_operators(SBPp, Ny, Nz, μ, Ly, Lz)
    # and factor main matrix with Cholesky
    M = cholesky(Symmetric(M̃))

    # initialize time and vector g that stores boundary data:
    t = 0
    g = zeros((Ny+1) * (Nz+1))

    # initialize slip:
    δ = zeros(Nz+1)
  
    # initialize boundary data on 4 corners of domain:
    # Dirichlet on left/right
    bc_Dirichlet = (lf, y, z) -> (2-lf) .* (0 * y .+ 0 .* z) + (lf-1) .* (0 .* y .+ 0 .* z)
    # traction free (Neumann) on top and bottom.  y = 0 is Earth's free surface.
    bc_Neumann   = (lf, y, z, ny, nz) -> zeros(size(y))
    # modify boundary data vector g with data from four sides:
    bdry_vec_mod!(g, F, τ, y, z, bc_Dirichlet, bc_Neumann)

    # Solve the linear system to get initial displacement vector u:
    u = M \ g
  
    # create the computational vector zf that is fault coordinates:
     zf = z   # i.e. just z-variable here
     (mm, δNp) = findmin(abs.(RSWf .- zf))  # calculate width of of RS fault 
     @assert zf[δNp] ≈ RSWf

    # initialize shear stress change:
    Δτ = zeros(Nz+1)

    # Assemble fault variables/data:
    RSa = zeros(δNp)
    for n = 1:δNp
        RSa[n] = RSamin - (RSamin - RSamax) *
          min(1, max(0, (RSH1 - zf[n])/(RSH1 - RSH2)))
    end

  
    # Calculate prestress:
    τz0 = σn * RSamax * asinh(RSVinit / (2 * RSV0) *
                                 exp.((RSf0 + RSb * log.(RSV0 / RSVinit)) /
                                      RSamax)) + η * RSVinit

    # initialize the state variable θ:
    θ = (RSDc ./ RSV0) .* exp.((RSa ./ RSb) .* log.((2 .* RSV0 ./ RSVinit) .*
      sinh.((τz0 .- η .* RSVinit) ./ (RSa .* σn))) .- RSf0 ./ RSb)

    # calculate scaled state variable ψ (for numerical purposes):
    ψ0 = RSf0 .+ RSb .* log.(RSV0 .* θ ./ RSDc)

    # initialize vector sent to ODE [ψ, δ]:
    ψδ = zeros(δNp + Nz + 1)  #because length(ψ) = δNp,  length(δ) = N+1
    # fill vector:  
    ψδ[1:δNp] .= ψ0
    ψδ[δNp+1:δNp + Nz + 1] .= δ


    function find_station_index(stations, grid_points)
      numstations = length(stations)
      station_ind = zeros(numstations)
      for i in range(1, stop=numstations)
        station_ind[i] = argmin(abs.(grid_points .- stations[i]))
      end
      return Integer.(station_ind)
    end

    # Stations to record on-fault time series:
    stations = [0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 25, 30, 35] # physical stations (km)
    station_indices = find_station_index(stations, zf)

    # Fault locations to record slip evolution output:
    flt_loc = 0:0.5:RSWf # physical stations every 500 m (km)
    flt_loc_indices = find_station_index(flt_loc, zf)

    u_old = copy(u)

  # set up parameters sent to right hand side of ODE:
  odeparam = (reject_step = [false],
              Vp=Vp,
              M = M,
              u=u,
              u_old=u_old,
              Δτ = Δτ,
              g = g,
              μshear=μshear,
              RSa=RSa,
              RSb=RSb,
              σn=σn,
              η=η,
              RSV0=RSV0,
              τz0=τz0,
              RSDc=RSDc,
              RSf0=RSf0,
              δNp = δNp,
              N = Nz,
              F = F,
              y = y ,
              z = z,
              τ = τ,
              HfI_FT = HfI_FT,
              save_stride_fields = stride_time # save every save_stride_fields time steps
             )

  #dψV = zeros(δNp + Ns + 1)
  tspan = (0, sim_years * year_seconds)

  # Set up the ODE problem to be solved, with rhs defined by odefun, 
  # on timespan tspan with initial condition ψδ and parameters odeparam:
  prob = ODEProblem(odefun, ψδ, tspan, odeparam)
  
  # this function gets called within rhs to enforce additional stability: 
  function stepcheck(_, p, _)
    if p.reject_step[1]
      p.reject_step[1] = false
      println("reject")
      return true
    end
    return false
  end

  #ODEresults = ODE_results([], [], [], Dict(i => [] for i = 1:length(stations)))

  # devol.txt is a file that stores time, max(V) and slip at all the stations:
  open("devol.txt", "w") do io
    write(io,"0.0 0.0 ")
        for i in 1:length(flt_loc)
          write(io,"$(flt_loc[i]) ")
        end
        write(io,"\n")
  end

  #write out initial data into devol.txt:
  vv = Array{Float64}(undef, 1, 2+length(flt_loc))
      vv[1] = t
      vv[2] = log10(RSVinit)
      vv[3:end] = δ[flt_loc_indices]
      open("devol.txt", "a") do io
         writedlm(io, vv)
      end


  # cb_mod gets called after every successful time step computed in ODE solver
  # (i.e. slip is written out to the text file devol.txt):
  cb_mod = SavingCallback((ψδ, t, i) -> write_text_slip(ψδ, t, i, zf, flt_loc, flt_loc_indices, odeparam, "BP1_", 10 * year_seconds), SavedValues(Float64, Float64))

  # Solve the ODE problem "prob" with Tsit5 (a Runge-Kutta method):
  sol = solve(prob, Tsit5(); dt=0.01,
              atol = 1e-14, rtol = 1e-14, save_everystep=true, gamma = 0.2,
              internalnorm=(x, _)->norm(x, Inf), callback=cb_mod)

 
  return (sol, zf, δNp)
end



# animate_slip will plot slip profiles against depth for every time step computed:
function animate_slip(S, δNp, zf, stride_time)

  m = length(zf)
  no_time_steps = size(S.t)
  slip_final = S.u[end][end]

  for i = 1:stride_time:no_time_steps[1]

    slip_t = S.u[i][δNp+1:end] # slip at time t
    #pyplot()
    display(plot(slip_t, -zf, xtickfont=font(18),
    ytickfont=font(18),
    guidefont=font(18),
    legendfont=font(18), ylabel = "Depth (km)", xlabel = "Slip (m)", xlims = (0, slip_final)))
    sleep(0.1)
  end

  #nothing
end


function interp1(xpt, ypt, x)

      knots = (xpt,) 
      itp = interpolate(knots, ypt, Gridded(Linear()))
      itp[x]  # endpoints of x must be between xpt[1] and xpt[end]
end


# find_ind() differentiates b/t phases by defining
# interseismic when max slip rate < 10^-3 m/s
# mv is maximum slip rate (log10 m/s) 
function find_ind(mv)
  ind = [1]
  int = 1
  cos = 0
  for i = 2:length(mv)
    if mv[i] > -3 && int == 1 && cos == 0
      append!(ind, i);
      int = 0;
      cos = 1;
    end
  
    if mv[i] < -3 && int == 0 && cos == 1
      append!(ind, i-1)
      int = 1
      cos = 0
    end
  end


  ind = append!(ind, length(mv));  #tack on for plotting any part of an incomplete coseismic/interseismic phase
  
  return ind
end

# plot_slip will plot slip contours from devol.txt - every 5 years in blue during interseismic, 
# every 1 second in red during coseismic
function plot_slip(filename)

  grid = readdlm(filename, Float64)
  sz = size(grid)
  flt_loc = grid[1,3:end]
  T = grid[2:sz[1],1]
  maxV = grid[2:end, 2]
  slip = grid[2:sz[1], 3:sz[2]]
  N = size(slip)[2]


  ind = find_ind(maxV);        #finds indices for inter/co-seismic phases
  interval = [5*31556926 1]   #plot every 5 years and every 1 second
  
  ct = 0   #this counts the number of events


  #Assumes an initial interseismic period
  #This for-loop only plots completed phases
  for i = 1:2:length(ind)-2
    
    T1 = T[ind[i]]:interval[1]:T[ind[i+1]];

    W1 = interp1(T,slip[:,1],T1)';
    
    for j = 2:N 
      w1 = interp1(T,slip[:,j],T1)';
      W1 = [W1; w1]
    end

    if i == 1
      plot(W1, -flt_loc, linecolor = :blue, legend = false) #interseismic phase
    else
      plot!(W1, -flt_loc, linecolor = :blue, legend = false) #interseismic phase
    end

   
    T1 = T[ind[i+1]]:interval[2]:T[ind[i+2]];


    W1 = interp1(T,slip[:,1],T1)';
    for j = 2:N 
      w1 = interp1(T,slip[:,j],T1)';
      W1 = [W1; w1]
    end

    plot!(W1, -flt_loc, linecolor = :red, legend = false) #interseismic phase

    ct = ct+1;
  end

  
  # plot remainder of an incomplete interseismic period:
  i = length(ind)-1;
  T1 = T[ind[i]]:interval[1]:T[ind[i+1]];
  W1 = interp1(T,slip[:,1],T1)';
      
      for j = 2:N 
        w1 = interp1(T,slip[:,j],T1)';
        W1 = [W1; w1]
      end

      plot!(W1, -flt_loc, linecolor = :blue, legend = false) #interseismic phase

end


# call main function: 
(S, zf, δNp) = main()

# animate slip (uncomment if desired):
# animate_slip(S, δNp, zf, 10)

# plot slip (uncomment if desired):
# plot_slip("devol.txt")