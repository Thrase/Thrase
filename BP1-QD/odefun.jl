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


function setupfaultstations(locations)
  T = eltype(locations)
  @assert size(locations,2) == 2
end


struct ODE_results
  t_list::Array{Any,1}
  V_list::Array{Any,1}
  δ_list::Array{Any,1}
 #  station_1::Array{Any,1}
 #  station_2::Array{Any,1}
 #  station_3::Array{Any,1}
 #  station_4::Array{Any,1}
 #  station_5::Array{Any,1}
 #  station_6::Array{Any,1}
 #  station_7::Array{Any,1}
 #  station_8::Array{Any,1}
 #  station_9::Array{Any,1}
 #  station_10::Array{Any,1}
 #  station_11::Array{Any,1}
 stations::Dict{Int64,Array{Any,1}}
end

function stopping_criteria(ψδ,t,i)
 Vmax = 0.0
 if isdefined(i,:fsallast)
   δNp = div(length(ψδ),2)
   dψV = i.fsallast
   dψ = @view dψV[1:δNp]
   V = @view dψV[δNp .+ (1:δNp)]
   Vmax = maximum(abs.(extrema(V)))
   δ = @view ψδ[δNp .+ (1:δNp)]
   ψ = @view ψδ[1:δNp]
   return Vmax >= 1e-3
 end
end

affect!(integrator)= terminate!(integrator)

function savestop(ψδ,t,i,ODEresults,p)
  Vmax = 0.0

  if isdefined(i,:fsallast)
  
    δNp = p.δNp
    N = p.N
    dψV = i.fsallast
    dψ = @view dψV[1:δNp]
    push!(ODEresults.t_list,t)

    
    V = @view dψV[δNp .+ (1:N+1)]
    Vmax = maximum(abs.(extrema(V)))
    δ = @view ψδ[δNp .+ (1:δNp)]
    ψ = @view ψδ[1:δNp]
    # p.u_old .= p.u
    if Vmax >= 1e-3
    
      push!(p.counter,1)
      if length(p.counter) >= 2
        terminate!(i)
        open("data_for_toby.txt","w") do io
          write(io,"$(ODEresults.t_list[end]) $(ODEresults.t_list[end-1])\n")
          write(io,"*"^40,"\n")
          writedlm(io,ψ,' ')
          write(io,"*"^40,"\n")
          writedlm(io,[p.u p.u_old],' ')
        end
      end
      p.u_old .= p.u
    end
 
 
  end

  Vmax
 end



function write_text_slip(ψδ,t,i,yf,stations,station_indices,p,base_name="",tdump=100)
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
    
    if mod(ctr[], p.save_stride_fields) == 0 || t == (sim_years ./ 31556926)
      vv = Array{Float64}(undef, 1, 2+length(stations))
      vv[1] = t
      vv[2] = log10(Vmax)
      vv[3:end] = δ[station_indices]
    
      open("devol.txt", "a") do io
        
        writedlm(io, vv)
      end
      
    end
    global ctr[] += 1
    @show ctr[]
  end

     Vmax
 end
