using Plots

function fault_params(fc, μ, cs)

    
    Wf = 40
    Hvw = 15
    Ht = 3
    δNp = findmin(abs.(Wf .- fc))[2]
    VWp = findmin(abs.((Hvw + Ht) .- fc))[2]
    a0 = .010
    amax = .025
    b = .015
    η = μ / (2*cs)
    
    function a_fun(y)
        if 0 <= y < Hvw
            return a0
        end
        if Hvw <= y < Hvw + Ht
            return a0 + (amax - a0)*(y-Hvw)/Ht
        end
        if Hvw + Ht <= y < Wf
            return amax
        end
        return amax
    end

    Dc = .008
    f0 = .6
    V0 = 1e-6
    Vp = 1e-9
    Vinit = 1e-9
    σn = 50.0
    a = a_fun.(fc[1:δNp])

    #plot(a, fc[1:δNp], yflip=true)
    #gui()
    #quit()
    
    τ⁰ = σn * amax * asinh(Vinit / (2 * V0) *
        exp.((f0 + b * log.(V0 / Vinit)) /
        amax))

    τq⁰ = τ⁰ + η * Vinit

        RS = (σn = σn,
              a = a,
              b = b,
              Dc = Dc,
              f0 = f0,
              V0 = V0,
              δNp = δNp,
              VWp = VWp,
              Vp = Vp,
              η = η,
              τ⁰ = τ⁰,
              τq⁰ = τq⁰)


    ψ0 = θtoψ((Dc ./ V0) .*
        exp.((a ./ b) .*
        log.((2 .* V0 ./ Vinit) .*
        sinh.((τ⁰) ./ (a .* σn))) .- f0 ./ b), RS)
    
    
    return RS, ψ0
    
end

ψtoθ(ψ, RS) = RS.Dc/RS.V0 .* exp.((ψ .- RS.f0)./RS.b)

θtoψ(θ,RS) = RS.f0 .+ RS.b * log.(RS.V0 * θ ./ RS.Dc)
