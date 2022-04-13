include("operators.jl")
include("material_params.jl")

using Printf

function Q_DYNAMIC_SOLVE!(dψV, ψδ, p, t)

    nn = p.nn
    δNp = p.RS.δNp
    Vp = p.RS.Vp
    Δτ = p.vars.Δτ
    τ = p.vars.τ
    M = p.ops.M̃
    u = p.vars.u
    ge = p.vars.ge
    RS = p.RS
    η = p.RS.η
    fc = p.metrics.facecoord[2][1]
    
    reject_step = p.reject_step
    if reject_step[1]
        return
    end
    
    ψ  = @view ψδ[(1:δNp)]
    δ =  @view ψδ[δNp .+ (1:nn)]
    dψ = @view dψV[(1:δNp)]
    V = @view dψV[δNp .+ (1:nn)]
    
    
    bc_Dirichlet = (lf, x, y) -> (2-lf)*(δ ./ 2) + (lf-1)*fill(t * Vp/2, size(x))
    bc_Neumann   = (lf, x, y) -> zeros(size(x))
    
    locbcarray_mod!(ge, p, bc_Dirichlet, bc_Neumann)
    
    u[:] = M \ ge
    

    Δτ .= - computetraction_mod(p, 1, u, δ)

    # solve for velocity point by point and set state derivative
    for n = 1:δNp
        ψn = ψ[n]
        an = RS.a[n]
        
        τn = Δτ[n] + RS.τq⁰
        τ[n] = τn
        
        if isnan(τn) || !isfinite(τn)
            #println("τ reject")
            reject_step[1] = true
            return
        end

        VR = abs(τn / η)
        VL = -VR
        Vn = V[n]
        obj_rs(V) = rateandstateQ(V, ψn, RS.σn, τn, RS.η, an, RS.V0)
        (Vn, _, iter) = newtbndv(obj_rs, VL, VR, Vn; ftol = 1e-12,
                                 atolx = 1e-12, rtolx = 1e-12)
        
        if isnan(Vn) || iter < 0 || !isfinite(Vn)
            #println("V reject")
            reject_step[1] = true
            return
        end

        V[n] = Vn
        dψ[n] = (RS.b * RS.V0 / RS.Dc) * (exp((RS.f0 - ψn) / RS.b) - abs(Vn) / RS.V0)

        if !isfinite(dψ[n]) || isnan(dψ[n])
            #println("ψ reject")
            dψ[n] = 0
            reject_step[1] = true
            return
        end
    end
    
    V[δNp+1:nn] .= RS.Vp

    nothing
end


# timestepping rejection function
function stepcheck(_, p, _)
    if p.reject_step[1]
        p.reject_step[1] = false
        return true
    end
    return false
end

# function for every accepted timstep with integrator stopping condition
function STOPFUN_Q(ψδ,t,i)
    
    if isdefined(i,:fsallast)

        δNp = i.p.RS.δNp
        nn = i.p.nn
        RS = i.p.RS
        τ = i.p.vars.τ
        t_prev = i.p.vars.t_prev
        year_seconds = i.p.vars.year_seconds
        u_prev = i.p.vars.u_prev
        u = i.p.vars.u
        fault_coord = i.p.metrics.facecoord[2][1]
        Lw = i.p.Lw
        io = i.p.io
        η = i.p.RS.η
        pf = i.p.io.pf
        
        dψV = i.fsallast
        ψ = ψδ[(1:δNp)]
        δ = ψδ[δNp .+ (1:nn)]
        V = @view dψV[δNp .+ (1:nn)]
        Vmax = maximum(abs.(V))
        

            
        θ = ψtoθ(ψ, RS)
        
        write_out(δ, V, τ, θ,
                  fault_coord,
                  Lw,
                  t,
                  io.station_names,
                  η)
        
        
        if pf[1] % 30 == 0

            
            plot!(δ[1:δNp], fault_coord[1:δNp], yflip = true, ylabel="Depth",
                  xlabel="Slip", linecolor=:blue, linewidth=.1,
                  legend=false)
            gui()
            
            write_out_ss(δ, V, τ, ψ, t,
                         io.slip_file,
                         io.stress_file,
                         io.slip_rate_file,
                         io.state_file,
                         η)
           
        end

        year_count = (t - t_prev[2])/year_seconds

        
        if Vmax >= 1e-2 && year_count > (t_prev[2] + 20)
            return true
        end
        
        pf[1] += 1
        u_prev .= u
        t_prev[1] = t
        
    end
    
    return false
        
end



function DYNAMIC_SOLVE!(dq, q, p, t)

    nn = p.nn
    RS = p.RS
    Λ = p.ops.Λ
    sJ = p.metrics.sJ[1]
    Z̃f = p.ops.Z̃f[1]
    L = p.ops.L[1]
    H = p.ops.H[1]
    P̃I = p.ops.P̃I
    JIHP = p.ops.JIHP
    nCnΓ1 = p.ops.nCnΓ1
    nBBCΓL1 = p.ops.nBBCΓL1
    vf = p.vars.Δτ
    τ̃f = p.vars.τ
    v̂_fric = p.vars.v̂_fric
    δNp = p.RS.δNp
    Vp = p.RS.Vp
    
    u = @view q[1:nn^2]
    v = @view q[nn^2 + 1 : 2nn^2]
    û1 = @view q[2nn^2 + 1 : 2nn^2 + nn]
    ψ = @view q[2nn^2 + 4nn + 1 : 2nn^2 + 4nn + δNp]
    dψ = @view dq[2nn^2 + 4nn + 1 : 2nn^2 + 4nn + δNp]
    dû1 = @view dq[2nn^2 + 1 : 2nn^2 + nn]
    dv =  @view dq[nn^2 + 1 : 2nn^2]
    
    dq .= Λ[1 : 2nn^2 + 4nn + δNp, 1 : 2nn^2 + 4nn + δNp] * q
    # get velocity on fault
    vf .= L * v
    # compute numerical traction on face 1
    τ̃f .= nBBCΓL1 * u + nCnΓ1 * û1

    # Root find for RS friction
    for n in 1:nn

        if n <= δNp
        
            v̂_root(v̂) = rateandstateD(v̂,
                                      Z̃f[n],
                                      vf[n],
                                      sJ[n],
                                      ψ[n],
                                      RS.a[n],
                                      τ̃f[n],
                                      RS.τ⁰,
                                      RS.σn,
                                      RS.V0)

            left = vf[n] + ((sJ[n] * RS.τ⁰) - τ̃f[n])/Z̃f[n]
            right = -left
        
            if left > right  
                tmp = left
                left = right
                right = tmp
            end
            
            (v̂n, _, _) = newtbndv(v̂_root, left, right, vf[n]; ftol = 1e-12,
                                  atolx = 1e-12, rtolx = 1e-12)
            
            if isnan(v̂n)
                println("Not bracketing root")
            end
            
            dû1[n] = v̂n
            dv[1 + (n-1)*nn] += H[n,n] * Z̃f[n] * v̂n
            dψ[n] = (RS.b .* RS.V0 ./ RS.Dc) .*
                (exp.((RS.f0 .- ψ[n]) ./ RS.b) .-
                abs.(2 .* v̂n) ./ RS.V0)
            
        else
            v̂_fric[n] = Vp/2
            dv[1 + (n-1)*nn] += H[n,n] * τ̃f[n]
        end
    end

    dv .= JIHP * dv
    
end




# bracketed newton method
function newtbndv(func, xL, xR, x; ftol = 1e-6, maxiter = 500, minchange=0,
                  atolx = 1e-4, rtolx = 1e-4)
    (fL, _) = func(xL)
    (fR, _) = func(xR)
    if fL .* fR > 0
        return (typeof(x)(NaN), typeof(x)(NaN), -maxiter)
    end

    (f, df) = func(x)
    dxlr = xR - xL

    for iter = 1:maxiter
        dx = -f / df
        x  = x + dx

        if x < xL || x > xR || abs(dx) / dxlr < minchange
            x = (xR + xL) / 2
            dx = (xR - xL) / 2
        end

        (f, df) = func(x)

        if f * fL > 0
            (fL, xL) = (f, x)
        else
            (fR, xR) = (f, x)
        end
        dxlr = xR - xL

        if abs(f) < ftol && abs(dx) < atolx + rtolx * (abs(dx) + abs(x))
            return (x, f, iter)
        end
    end
    return (x, f, -maxiter)
end


function D_timestep!(q, f!, p, dt, (t0, t1), Δq = similar(q), Δq2 = similar(q))
    
    T = eltype(q)
    nn = p.nn
    L = p.ops.L[1]
    VWp = p.RS.VWp
    δNp = p.RS.δNp
    
    pf = p.io.pf
    RS = p.RS
    
    τ̃ = p.vars.τ[1:δNp]
    sJ = p.metrics.sJ[1][1:δNp]
    Z̃f = p.ops.Z̃f[1][1:δNp]
    fc = p.metrics.facecoord[2][1][1:δNp]
    Lw = p.Lw
    io = p.io
    pf = p.io.pf
    d_to_s = p.d_to_s
    
    v = @view q[nn^2 + 1:2nn^2]
    u = @view q[1:nn^2]
    ψ = @view q[2nn^2 + 4*nn + 1 : 2nn^2 + 4*nn + δNp]
    v̂ = p.vars.ge
    
    RKA = (
        T(0),
        T(-567301805773 // 1357537059087),
        T(-2404267990393 // 2016746695238),
        T(-3550918686646 // 2091501179385),
        T(-1275806237668 // 842570457699),
    )

    RKB = (
        T(1432997174477 // 9575080441755),
        T(5161836677717 // 13612068292357),
        T(1720146321549 // 2090206949498),
        T(3134564353537 // 4481467310338),
        T(2277821191437 // 14882151754819),
    )

    RKC = (
        T(0),
        T(1432997174477 // 9575080441755),
        T(2526269341429 // 6820363962896),
        T(2006345519317 // 3224310063776),
        T(2802321613138 // 2924317926251),
    )

    nstep = ceil(Int, (t1 - t0) / dt)
    dt = (t1 - t0) / nstep

    fill!(Δq, 0)
    fill!(Δq2, 0)

    pf[1] = .01
    pf[2] = .5
    
    for step in 1:nstep
        t = t0 + (step - 1) * dt
        if mod(nstep, 100) == 0
        end
        
        for s in 1:length(RKA)
            f!(Δq2, q, p, t + RKC[s] * dt)
            v̂[1:δNp] .= Δq2[2nn^2 + 1 : 2nn^2 + δNp]
            Δq .+= Δq2
            q .+= RKB[s] * dt * Δq
            Δq .*= RKA[s % length(RKA) + 1]
        end
        
        
 
        θ = ψtoθ(ψ,RS)
        vf = (L*v)[1:δNp]
        τ̂ = -τ̃ ./ sJ .- Z̃f .* (v̂[1:δNp] .- vf) ./ sJ .+ RS.τ⁰

        if step == ceil(Int, pf[1]/dt)
            
            write_out(2*(L*u)[1:δNp],
                      2*v̂[1:δNp], τ̂, θ,
                      fc,
                      Lw, t,
                      io.station_names)
            
            pf[1] += .01
        end
        
        if step == ceil(Int, pf[2]/dt)

            write_out_ss(2*(L*u)[1:δNp],
                         2v̂[1:δNp], ψ, τ̂, t,
                         io.slip_file,
                         io.stress_file,
                         io.slip_rate_file,
                         io.state_file)

            plot!(2*(L*u)[1:δNp], fc, yflip = true, ylabel="Depth",
                  xlabel="Slip", linecolor=:red, linewidth=.1,
                  legend=false)
            gui()
            
            pf[2] += .5
        end

        if (2 * maximum(v̂[1:VWp])) < d_to_s && pf[1] > 7.0
            return t
        end

    end
    return t0
end
