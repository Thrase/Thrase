using DelimitedFiles
using LinearAlgebra
using SparseArrays
using Printf
include("DiagonalSBP.jl")

⊗(A,B) = kron(A,B)

function operators_dynamic(p, Nr, Ns, μ, ρ, R, faces, metrics, 
                     τscale = 2,
                     crr = metrics.crr,
                     css = metrics.css,
                     crs = metrics.crs)


    
    csr = crs
    J = metrics.J

    hr = 2/Nr
    hs = 2/Ns

    hmin = min(hr, hs)
    
    r = -1:hr:1
    s = -1:hs:1
    
    Nrp = Nr + 1
    Nsp = Ns + 1
    Np = Nrp * Nsp
    Nn = Np
    nn = Nrp
    # Derivative operators for the rest of the computation
    (Dr, HrI, Hr, r) = D1(p, Nr; xc = (-1,1))
    Qr = Hr * Dr
    QrT = sparse(transpose(Qr))

    (Ds, HsI, Hs, s) = D1(p, Ns; xc = (-1,1))
    Qs = Hs * Ds
    QsT = sparse(transpose(Qs))

    # Identity matrices for the comuptation
    Ir = sparse(I, Nrp, Nrp)
    Is = sparse(I, Nsp, Nsp)

    (_, S0e, SNe, _, _, Ae, _) = variable_D2(p, Nr, rand(Nrp))
    IArr = Array{Int64,1}(undef,Nsp * length(Ae.nzval))
    JArr = Array{Int64,1}(undef,Nsp * length(Ae.nzval))
    VArr = Array{Float64,1}(undef,Nsp * length(Ae.nzval))
    stArr = 0

    for j = 1:Nsp
        rng = (j-1) * Nrp .+ (1:Nrp)
        (_, S0e, SNe, _, _, Ae, _) =  variable_D2(p, Nr, crr[rng])
        (Ie, Je, Ve) = findnz(Ae)
        IArr[stArr .+ (1:length(Ve))] = Ie .+ (j-1) * Nrp
        JArr[stArr .+ (1:length(Ve))] = Je .+ (j-1) * Nrp
        VArr[stArr .+ (1:length(Ve))] = Hs[j,j] * Ve
        stArr += length(Ve)
    end

    Ãrr = sparse(IArr[1:stArr], JArr[1:stArr], VArr[1:stArr], Np, Np)

    (_, S0e, SNe, _, _, Ae, _) =  variable_D2(p, Ns, rand(Nsp))
    IAss = Array{Int64,1}(undef,Nrp * length(Ae.nzval))
    JAss = Array{Int64,1}(undef,Nrp * length(Ae.nzval))
    VAss = Array{Float64,1}(undef,Nrp * length(Ae.nzval))
    stAss = 0

    for i = 1:Nrp
        rng = i .+ Nrp * (0:Ns)
        (_, S0e, SNe, _, _, Ae, _) =  variable_D2(p, Ns, css[rng])

        (Ie, Je, Ve) = findnz(Ae)
        IAss[stAss .+ (1:length(Ve))] = i .+ Nrp * (Ie .- 1)
        JAss[stAss .+ (1:length(Ve))] = i .+ Nrp * (Je .- 1)
        VAss[stAss .+ (1:length(Ve))] = Hr[i,i] * Ve
        stAss += length(Ve)

    end

    Ãss = sparse(IAss[1:stAss], JAss[1:stAss], VAss[1:stAss], Np, Np)

    Ãsr = (QsT ⊗ Ir) * sparse(1:length(crs), 1:length(crs), view(crs, :)) * (Is ⊗ Qr)
    Ãrs = (Is ⊗ QrT) * sparse(1:length(csr), 1:length(csr), view(csr, :)) * (Qs ⊗ Ir)

    Ã = Ãrr + Ãss + Ãrs + Ãsr
   
    # volume quadrature
    H̃ = kron(Hr, Hs)
    H̃inv = spdiagm(0 => 1 ./ diag(H̃))
    
    # diagonal rho
    P̃ = spdiagm(0 => repeat([ρ], Nrp * Nsp))
    P̃inv = spdiagm(0 => (1 ./ repeat([ρ], Nrp * Nsp)))
    JI = spdiagm(0 => reshape(metrics.JI, Nrp*Nsp))
    
    er0 = sparse([1  ], [1], [1], Nrp, 1)
    erN = sparse([Nrp], [1], [1], Nrp, 1)
    es0 = sparse([1  ], [1], [1], Nsp, 1)
    esN = sparse([Nsp], [1], [1], Nsp, 1)

    er0T = sparse([1], [1  ], [1], 1, Nrp)
    erNT = sparse([1], [Nrp], [1], 1, Nrp)
    es0T = sparse([1], [1  ], [1], 1, Nsp)
    esNT = sparse([1], [Nsp], [1], 1, Nsp)

    cmax = maximum([maximum(crr), maximum(crs), maximum(css)])

    # Surface quadtrature matrices
    H1 = H2 = Hs 
    H3 = H4 = Hr

    H = (Hs, Hs, Hr, Hr)
    HI = (HsI, HsI, HrI, HrI)
    # Volume to Face Operators (transpose of these is face to volume)
    L = (convert(SparseMatrixCSC{Float64, Int64}, kron(Ir, es0)'),
         convert(SparseMatrixCSC{Float64, Int64}, kron(Ir, esN)'),
         convert(SparseMatrixCSC{Float64, Int64}, kron(er0, Is)'),
         convert(SparseMatrixCSC{Float64, Int64}, kron(erN, Is)'))
    # coefficent matrices
    Crr1 = spdiagm(0 => crr[1, :])
    Crs1 = spdiagm(0 => crs[1, :])
    Csr1 = spdiagm(0 => crs[1, :])
    Css1 = spdiagm(0 => css[1, :])
    
    Crr2 = spdiagm(0 => crr[Nrp, :])
    Crs2 = spdiagm(0 => crs[Nrp, :])
    Csr2 = spdiagm(0 => crs[Nrp, :])
    Css2 = spdiagm(0 => css[Nrp, :])
    
    Css3 = spdiagm(0 => css[:, 1])
    Crs3 = spdiagm(0 => crs[:, 1])
    Csr3 = spdiagm(0 => crs[:, 1])
    Crr3 = spdiagm(0 => crr[:, 1])
    
    Css4 = spdiagm(0 => css[:, Nsp])
    Crs4 = spdiagm(0 => crs[:, Nsp])
    Csr4 = spdiagm(0 => crs[:, Nsp])
    Crr4 = spdiagm(0 => crr[:, Nsp])

    (_, S0, SN, _, _) = D2(p, Nr, xc=(-1,1))[1:5]
    S0 = sparse(Array(S0[1,:])')
    SN = sparse(Array(SN[end, :])')
    
    # Boundars Derivatives
    B1r =  Crr1 * kron(Is, S0)
    B1s = Crs1 * L[1] * kron(Ds, Ir)
    B2r = Crr2 * kron(Is, SN)
    B2s = Crs2 * L[2] * kron(Ds, Ir)
    B3s = Css3 * kron(S0, Ir)
    B3r = Csr3 * L[3] * kron(Is, Dr)
    B4s = Css4 * kron(SN, Ir)
    B4r = Csr4 * L[4] * kron(Is, Dr)
    

    (xf1, xf2, xf3, xf4) = metrics.facecoord[1]
    (yf1, yf2, yf3, yf4) = metrics.facecoord[2]

    Z̃f = (metrics.sJ[1] .* sqrt(μ * ρ),
          metrics.sJ[2] .* sqrt.(ρ * μ),
          metrics.sJ[3] .* sqrt.(ρ * μ),
          metrics.sJ[4] .* sqrt.(ρ * μ))

    # Penalty terms
        if p == 2
        l = 2
        β = 0.363636363
        α = 1 / 2
        θ_R = 1.0
    elseif p == 4
        l = 4
        β = 0.2505765857
        α = 17 / 48
        θ_R = 0.5776
    elseif p == 6
        l = 7
        β = 0.1878687080
        α = 13649 / 43200
        θ_R = 0.3697
    else
        error("unknown order")
    end

    ψmin_r = reshape(crr, Nrp, Nsp)
    ψmin_s = reshape(css, Nrp, Nsp)
    @assert minimum(ψmin_r) > 0
    @assert minimum(ψmin_s) > 0
    
    hr = 2 / Nr
    hs = 2 / Ns

    ψ1 = ψmin_r[  1, :]
    ψ2 = ψmin_r[Nrp, :]
    ψ3 = ψmin_s[:,   1]
    ψ4 = ψmin_s[:, Nsp]
    
    for k = 2:l
        ψ1 = min.(ψ1, ψmin_r[k, :])
        ψ2 = min.(ψ2, ψmin_r[Nrp+1-k, :])
        ψ3 = min.(ψ3, ψmin_s[:, k])
        ψ4 = min.(ψ4, ψmin_s[:, Nsp+1-k])
    end
    
    τs1 = (2τscale / hr) * (crr[  1, :].^2 / β + crs[  1, :].^2 / α) ./ ψ1
    τs2 = (2τscale / hr) * (crr[Nrp, :].^2 / β + crs[Nrp, :].^2 / α) ./ ψ2
    τs3 = (2τscale / hs) * (css[:,   1].^2 / β + crs[:,   1].^2 / α) ./ ψ3
    τs4 = (2τscale / hs) * (css[:, Nsp].^2 / β + crs[:, Nsp].^2 / α) ./ ψ4
    
    
    τs = (sparse(1:Nsp, 1:Nsp, τs1),
          sparse(1:Nsp, 1:Nsp, τs2),
          sparse(1:Nrp, 1:Nrp, τs3),
          sparse(1:Nrp, 1:Nrp, τs4))

    τR1 = (1/(θ_R*hr))*Is
    τR2 = (1/(θ_R*hr))*Is
    τR3 = (1/(θ_R*hs))*Ir
    τR4 = (1/(θ_R*hs))*Ir
    
    p1 = ((crr[  1, :]) ./ ψ1)
    p2 = ((crr[Nrp, :]) ./ ψ2)
    p3 = ((css[:,   1]) ./ ψ3)
    p4 = ((css[:, Nsp]) ./ ψ4)
   
    P1 = sparse(1:Nsp, 1:Nsp, p1)
    P2 = sparse(1:Nsp, 1:Nsp, p2)
    P3 = sparse(1:Nrp, 1:Nrp, p3)
    P4 = sparse(1:Nrp, 1:Nrp, p4)

    # dynamic penalty matrices
    Γ = ((2/(α*hr))*Is + τR1 * P1,
         (2/(α*hr))*Is + τR2 * P2,
         (2/(α*hs))*Ir + τR3 * P3,
         (2/(α*hs))*Ir + τR4 * P4)

    JH = sparse(1:Np, 1:Np, view(J, :)) * (Hs ⊗ Hr)
    
    JIHP = JI * H̃inv * P̃inv
    
    Cf = ((Crr1, Crs1), (Crr2, Crs2), (Css3, Csr3), (Css4, Csr4))
    B = ((B1r, B1s), (B2r, B2s), (B3s, B3r), (B4s, B4r))
    nl = (-1, 1, -1, 1)
    

    F = ((nl[1] * H[1] * (B[1][1] + B[1][2]))' - L[1]'*H[1]*τs[1],
         (nl[2] * H[2] * (B[2][1] + B[2][2]))' - L[2]'*H[2]*τs[2],
         (nl[3] * H[3] * (B[3][1] + B[3][2]))' - L[3]'*H[3]*τs[3],
         (nl[4] * H[4] * (B[4][1] + B[4][2]))' - L[4]'*H[4]*τs[4])

    HIFT = ((nl[1] * (B[1][1] + B[1][2])) - (τs[1] ⊗ er0'),
            (nl[2] * (B[2][1] + B[2][2])) - (τs[2] ⊗ er0'),
            (nl[3] * (B[3][1] + B[3][2])) - (τs[3] ⊗ er0'),
            (nl[4] * (B[4][1] + B[4][2])) - (τs[4] ⊗ er0'))
    
    M̃ = copy(Ã)
        #M̃ = Ã
        for i in 1:4
            M̃ .+= -L[i]' * (nl[i] * H[i] * (B[i][1] + B[i][2])) -
                (nl[i] * H[i] * (B[i][1] + B[i][2]))' * L[i] +
                L[i]'*H[i]*τs[i]*L[i]
        
            if i == 3 || i == 4
                M̃ -= F[i] * (Diagonal(1 ./ (diag(τs[i]))) * HI[i]) * F[i]'
            end
        end  
    
    # accleration blocks
    dv_u = -Ã
    
    for i in 1:4
        if faces[i] == 0
            dv_u .+= (L[i]' * H[i] * (nl[i] * (B[i][1] + B[i][2]) - Cf[i][1] * Γ[i] * L[i])) +
                nl[i] * (B[i][1]' + B[i][2]') * H[i] * L[i]
        else
            dv_u .+= (L[i]' * H[i] * ((1 - R[i])/2 .* (nl[i] * (B[i][1] + B[i][2]) - Cf[i][1] * Γ[i] * L[i]))) +
                nl[i] * (B[i][1]' + B[i][2]') * H[i] * L[i]
        end
    end

    dv_v = spzeros(Nn, Nn)
    for i in 1:4
        if faces[i] == 0
            dv_v .+=  L[i]' * H[i] * (-Z̃f[i] .* L[i])
        else
            dv_v .+=  L[i]' * H[i] * (-(1 - R[i])/2 .* Z̃f[i] .* L[i])
        end
    end

    dv_û = spzeros(Nn, 4nn)
    dû_u = spzeros(4nn, Nn)
    dû_v = spzeros(4nn, Nn)

    for i in 1:4
        if faces[i] == 0
                        dv_û[ : , (i-1) * nn + 1 : i * nn] .=
                (L[i]' * H[i] * Cf[i][1] * Γ[i]) -
                nl[i] * (B[i][1]' + B[i][2]') * H[i]
        else
            dv_û[ : , (i-1) * nn + 1 : i * nn] .=
                (L[i]' * H[i] * ((1 - R[i])/2 .* Cf[i][1] * Γ[i])) -
                nl[i] * (B[i][1]' + B[i][2]') * H[i]
        
            dû_u[(i-1) * nn + 1 : i * nn , : ] .=
                -(1 + R[i])/2 .* (nl[i] * (B[i][1] + B[i][2]) -
                Cf[i][1] * Γ[i] * L[i])./Z̃f[i]
        
            dû_v[(i-1) * nn + 1 : i * nn , : ] .= (1 + R[i])/2 .* L[i]
        end
    end

    
    dû_û = spzeros(4nn, 4nn)
    for i in 1:4
        if faces[i] != 0
            dû_û[(i-1) * nn + 1 : i * nn, (i-1) * nn + 1 : i * nn] .= -(1 + R[i])/2 .* (Cf[i][1] * Γ[i])./Z̃f[i]
        end
    end

    dû_ψ = spzeros(4nn, nn)

    Λ = [ spzeros(Nn, Nn) sparse(I, Nn, Nn) spzeros(Nn, 5nn)
          dv_u dv_v dv_û spzeros(Nn, nn)
          dû_u dû_v dû_û dû_ψ
          spzeros(nn, 2Nn + 5nn) ]

    nCnΓ1 = Crr1 * Γ[1]
    nBBCΓL1 = nl[1] * (B[1][1] + B[1][2]) - nCnΓ1 * L[1]
    
    (Λ = Λ,
     M̃ = cholesky(Symmetric(M̃)),
     F = F,
     HIFT = HIFT,
     τQ = (diag(τs[1]),
           diag(τs[2]),
           diag(τs[3]),
           diag(τs[4])),
     P̃I = P̃inv,
     H̃I = H̃inv,
     JI = JI,
     JIHP = JIHP,
     nCnΓ1 = nCnΓ1,
     nBBCΓL1 = nBBCΓL1,
     hmin = hmin,
     cmax = cmax,
     JH = JH,
     sJ = metrics.sJ,
     nx = metrics.nx,
     ny = metrics.ny,
     L = L,
     H = H,
     Z̃f = Z̃f)
    
end     


function locbcarray_mod!(ge, p, bc_Dirichlet, bc_Neumann)
    F = p.ops.F
    (xf, yf) = p.metrics.facecoord
    sJ = p.metrics.sJ
    τ = p.ops.τQ
    ge[:] .= 0
    
    for i in 1:4
        if i == 1 || i == 2
            vf = bc_Dirichlet(i, xf[i], yf[i])
        elseif i == 3 || i == 4
            gN = bc_Neumann(i, xf[i], yf[i])
            vf = sJ[i] .* gN ./ τ[i]
        end
        ge[:] -= F[i] * vf
    end
    
end

function computetraction_mod(p, f, u, δ)
    HIFT = p.ops.HIFT[f]
    τf = p.ops.τQ[f]
    sJ = p.metrics.sJ[f]
    return (HIFT * u + τf .* (δ .- δ / 2)) ./ sJ
end
    
function rateandstateQ(V, psi, σn, τ, η, a, V0)
    Y = (1 ./ (2 .* V0)) .* exp.(psi ./ a)
    f = a .* asinh.(V .* Y)
    dfdV  = a .* (1 ./ sqrt.(1 + (V .* Y).^2)) .* Y
    g    = σn .* f    + η .* V - τ
    dgdV = σn .* dfdV + η
    (g, dgdV)
end

function rateandstateD(v̂, z̃, v, sJ, ψ, a, τ̃, τ⁰, σn, V0)
    
    Y = (1 / (2 * V0)) * exp(ψ / a)
    f = a * asinh.(2v̂ * Y)
    dfdv̂  = a * (1 / sqrt(1 + (2v̂ * Y)^2)) * 2Y
    g = sJ * σn * f - sJ * τ⁰ + τ̃ + z̃*(v̂ - v)
    dgdv̂ = z̃ + sJ * σn * dfdv̂
    
    return (g, dgdv̂)
end


