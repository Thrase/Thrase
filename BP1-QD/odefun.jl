const year_seconds = 31556926
const sim_years = 1500
global const ctr = Ref{Int64}(1) 

using OrdinaryDiffEq
using DiffEqCallbacks
using Printf
using Plots
using DelimitedFiles

function odefun(dψV, ψδ, p, t)
  
  # unpack structure of parameters:
  Vp = p.Vp
  M = p.M
  u = p.u
  Δτ = p.Δτ
  g = p.g
  μshear = p.μshear
  RSa = p.RSa
  RSb = p.RSb
  σn = p.σn
  η = p.η
  RSV0 = p.RSV0
  τz0 = p.τz0
  RSDc = p.RSDc
  RSf0 = p.RSf0
  δNp = p.δNp
  N = p.N
  F = p.F
  τf = p.τf
  τ = p.τ
  y = p.y 
  z = p.z
  HfI_FT = p.HfI_FT


  # show current time:
  @show t ./ 31556926

  # unpack vector into ψ and δ:
  ψ  = @view ψδ[        (1:δNp) ]
  δ  = ψδ[ δNp .+ (1:N+1) ]

  # update functions for calculating boundary data:
  bc_Dirichlet = (lf, y, z) -> (2-lf)*(δ ./ 2) + (lf-1)*fill(t .* Vp./2, size(z))
  bc_Neumann   = (lf, y, z, ny, nz) -> zeros(size(y))
  bdry_vec_mod!(g, F, τ, y, z, bc_Dirichlet, bc_Neumann)

  # solve for displacements everywhere in domain
  u[:] = M \ g

  # set up rates of change for  state and slip and intialize to 0:
  dψ = @view dψV[       (1:δNp) ]
  V  = @view dψV[ δNp .+ (1:N+1)]
  dψ .= 0
  V  .= 0


  # Update the fault data
  Δτ .= 0
  lf1 = 1 # recall that face 1 (lf1) is face 1

  # compute change in shear stress:
  Δτ .= -computetraction_stripped(HfI_FT, τ, lf1, u, δ)

  # Loop over fault indices and update rates of change:
  for n = 1:δNp
    ψn = ψ[n]
    an = RSa[n]

    τn = Δτ[n] + τz0
    τf[n] = τn # store for writing out to file.
  

    VR = abs(τn / η)
    VL = -VR
    Vn = V[n]
    obj_rs(V) = rateandstate(V, ψn, σn, τn, η, an, RSV0)
    (Vn, _, iter) = newtbndv(obj_rs, VL, VR, Vn; ftol = 1e-9,
                                 atolx = 1e-9, rtolx = 1e-9)


    
    V[n] = Vn

    dψ[n] = (RSb * RSV0 / RSDc) * (exp((RSf0 - ψn) / RSb) - abs(Vn) / RSV0)
   
  end

  V[δNp+1:N+1] .= Vp


  nothing
end




function write_to_file(ψδ, t, i, zf,flt_loc, flt_loc_indices, stations, station_indices, p, base_name="", tdump=100)
  Vmax = 0.0

  if isdefined(i,:fsallast) 
    δNp = p.δNp
    N = p.N
    dψV = i.fsallast
    dψ = @view dψV[1:δNp]
    V = @view dψV[δNp .+ (1:N+1)]
    Vmax = maximum(abs.(extrema(V)))
    δ = @view ψδ[δNp .+ (1:δNp)]
    ψ = @view ψδ[1:δNp]
    τf = p.τf
 
    θ = (p.RSDc * exp.((ψ .- p.RSf0) ./ p.RSb)) / p.RSV0  # Invert ψ for θ.
    
    if mod(ctr[], p.save_stride_fields) == 0 || t == (sim_years ./ 31556926)
      vv = Array{Float64}(undef, 1, 2+length(flt_loc))
      vv[1] = t
      vv[2] = log10(Vmax)
      vv[3:end] = δ[flt_loc_indices]
      open("devol.txt", "a") do io
        writedlm(io, vv)
      end

      stations = ["000", "025", "050", "075", "100", "125", "150", "175", "200", "250", "300", "350"]
     
      for i = 1:length(station_indices)
        ww = Array{Float64}(undef, 1, 5)
        ww[1] = t
        ww[2] = δ[station_indices[i]]
        ww[3] = log10(V[station_indices[i]])
        ww[4] =  τf[station_indices[i]]
        ww[5] = log10(θ[station_indices[i]])
        open("fltst_dp"*stations[i]*".txt", "a") do io
            writedlm(io, ww)
        end
      end
    end
  
      
  
    global ctr[] += 1
    @show ctr[]
  

  end
     Vmax

  
  
end

 

      
function create_text_files(flt_loc, flt_loc_indices, stations, station_indices, t, RSVinit, δ, τz0, θ)

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

  # write out initial data into station files:

  # fltst_dpXXX.txt is a file that stores time and time-series of slip, log10(slip_rate), 
  # shear_stress and log10(state) at depth of z = XXX km, where XXX is each of the fault station depths.
  # First we write out initial data into each fltst_dpXXX.txt:

  stations = ["000", "025", "050", "075", "100", "125", "150", "175", "200", "250", "300", "350"]
  for n = 1:length(stations)
    XXX = "fltst_dp"*stations[n]*".txt"
    ww = Array{Float64}(undef, 1, 5)
    ww[1] = t
    ww[2] = δ[station_indices[n]]
    ww[3] = log10(RSVinit)
    ww[4] =  τz0
    ww[5] = log10(θ[station_indices[n]])
    open(XXX, "w") do io
        writedlm(io, ww)
    end
  end

end
