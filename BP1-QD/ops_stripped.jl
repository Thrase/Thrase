include("diagonal_sbp.jl")
using SparseArrays
using LinearAlgebra

function get_operators(p, Nr, Ns, μ, Ly, Lz)

    # In this project, take r = y, s = z (i.e. problem is set up for no coordinate transformation)
    Nrp = Nr + 1  
    Nsp = Ns + 1

    # "coefficient" matrices 
    crr = μ * ones(Nrp, Nsp)
    css = μ * ones(Nrp, Nsp)
    crs = zeros(Nrp, Nsp)
    csr = crs

# Derivative operators for the rest of the computation
    (Dr, HrI, Hr, r) = diagonal_sbp_D1(p, Nr; xc = (0,Ly))
    Qr = Hr * Dr
    QrT = sparse(transpose(Qr))

    (Ds, HsI, Hs, s) = diagonal_sbp_D1(p, Ns; xc = (0,Lz))
    Qs = Hs * Ds
    QsT = sparse(transpose(Qs))

    ⊗(A,B) = kron(A, B)
    Ir = sparse(I, Nrp, Nrp)
    Is = sparse(I, Nsp, Nsp)
    (_, S0r, SNr, _, _, Ar, _) = variable_diagonal_sbp_D2(p, Nr, μ * ones(Nrp); xc = (0,Ly))
    (_, S0s, SNs, _, _, As, _) = variable_diagonal_sbp_D2(p, Ns, μ * ones(Nsp); xc = (0,Lz))
  

    Er0 = sparse([1], [1], [1], Nrp, Nrp)
    ErN = sparse([Nrp], [Nrp], [1], Nrp, Nrp)
    Es0 = sparse([1], [1], [1], Nsp, Nsp)
    EsN = sparse([Nsp], [Nsp], [1], Nsp, Nsp)

    er0 = sparse([1  ], [1], [1], Nrp, 1)
    erN = sparse([Nrp], [1], [1], Nrp, 1)
    es0 = sparse([1  ], [1], [1], Nsp, 1)
    esN = sparse([Nsp], [1], [1], Nsp, 1)

    er0T = sparse([1], [1  ], [1], 1, Nrp)
    erNT = sparse([1], [Nrp], [1], 1, Nrp)
    es0T = sparse([1], [1  ], [1], 1, Nsp)
    esNT = sparse([1], [Nsp], [1], 1, Nsp)

        #
    # Store coefficient matrices as matrices
    #
    crs0 = sparse(Diagonal(crs[1:Nrp]))
    crsN = sparse(Diagonal(crs[Nrp*Ns .+ (1:Nrp)]))
    csr0 = sparse(Diagonal(csr[1   .+ Nrp*(0:Ns)]))
    csrN = sparse(Diagonal(csr[Nrp .+ Nrp*(0:Ns)]))

    #
    # Surface mass matrices
    #
    H1 = Hs
    H1I = HsI

    H2 = Hs
    H2I = HsI

    H3 = Hr
    H3I = HrI

    H4 = Hr
    H4I = HrI

    #
    # Penalty terms
    #
    if p == 2
        l = 2
        β = 0.363636363
        α = 1 / 2
    elseif p == 4
        l = 4
        β = 0.2505765857
        α = 17 / 48
    elseif p == 6
        l = 7
        β = 0.1878687080
        α = 13649 / 43200
    else
        error("unknown order")
    end

    ψmin = reshape((crr + css - sqrt.((crr - css).^2 + 4crs.^2)) / 2, Nrp, Nsp)
    @assert minimum(ψmin) > 0

    hr = Ly / Nr
    hs = Lz / Ns

    ψ1 = ψmin[  1, :]
    ψ2 = ψmin[Nrp, :]
    ψ3 = ψmin[:,   1]
    ψ4 = ψmin[:, Nsp]
    for k = 2:l
        ψ1 = min.(ψ1, ψmin[k, :])
        ψ2 = min.(ψ2, ψmin[Nrp+1-k, :])
        ψ3 = min.(ψ3, ψmin[:, k])
        ψ4 = min.(ψ4, ψmin[:, Nsp+1-k])
    end
    τscale = 2

    τ1 = (2τscale / hr) * (crr[  1, :].^2 / β + crs[  1, :].^2 / α) ./ ψ1
    τ2 = (2τscale / hr) * (crr[Nrp, :].^2 / β + crs[Nrp, :].^2 / α) ./ ψ2
    τ3 = (2τscale / hs) * (css[:,   1].^2 / β + crs[:,   1].^2 / α) ./ ψ3
    τ4 = (2τscale / hs) * (css[:, Nsp].^2 / β + crs[:, Nsp].^2 / α) ./ ψ4

    τ1 = sparse(1:Nsp, 1:Nsp, τ1)
    τ2 = sparse(1:Nsp, 1:Nsp, τ2)
    τ3 = sparse(1:Nrp, 1:Nrp, τ3)
    τ4 = sparse(1:Nrp, 1:Nrp, τ4)

    Ãrr = Hs ⊗ Ar
    Ãss = As ⊗ Hr
    Ã = Ãrr + Ãss 

    Sr0 = Hs ⊗ S0r
    SrN = Hs ⊗ SNr
    Ss0 = S0s ⊗ Hr
    SsN = SNs ⊗ Hr

    Sr0T = Hs ⊗ sparse(transpose(S0r))
    SrNT = Hs ⊗ sparse(transpose(SNr))
    Ss0T = sparse(transpose(S0s)) ⊗ Hr
    SsNT = sparse(transpose(SNs)) ⊗ Hr

    C̃1 =  (Sr0 + Sr0T) + ((csr0 * Qs + QsT * csr0) ⊗ Er0) + ((τ1 * H1) ⊗ Er0)
    C̃2 = -(SrN + SrNT) - ((csrN * Qs + QsT * csrN) ⊗ ErN) + ((τ2 * H2) ⊗ ErN)
    C̃3 =  (Ss0 + Ss0T) + (Es0 ⊗ (crs0 * Qr + QrT * crs0)) + (Es0 ⊗ (τ3 * H3))
    C̃4 = -(SsN + SsNT) - (EsN ⊗ (crsN * Qr + QrT * crsN)) + (EsN ⊗ (τ4 * H4))

    # TODO: Fix minus sign (reverse of the paper)
    G1 = -(Is ⊗ er0T) * Sr0 - ((csr0 * Qs) ⊗ er0T)
    G2 = +(Is ⊗ erNT) * SrN + ((csrN * Qs) ⊗ erNT)
    G3 = -(es0T ⊗ Ir) * Ss0 - (es0T ⊗ (crs0 * Qr))
    G4 = +(esNT ⊗ Ir) * SsN + (esNT ⊗ (crsN * Qr))

    F1 = G1' - ((τ1 * H1) ⊗ er0)
    F2 = G2' - ((τ2 * H2) ⊗ erN)
    F3 = G3' - (es0 ⊗ (τ3 * H3))
    F4 = G4' - (esN ⊗ (τ4 * H4))


    HfI_F1T = H1I * G1 - (τ1 ⊗ er0')
    HfI_F2T = H2I * G2 - (τ2 ⊗ erN')
    HfI_F3T = H3I * G3 - (es0' ⊗ τ3)
    HfI_F4T = H4I * G4 - (esN' ⊗ τ4)


    M̃ = Ã + C̃1 + C̃2 + C̃3 + C̃4

    H̃ = Hs ⊗ Hr
    # Modify the operator to handle the boundary conditions
    F = (F1, F2, F3, F4)
    τ = (τ1, τ2, τ3, τ4)
    HfI = (H1I, H2I, H3I, H4I)

    # Modify operators for the BC (Neumann only only faces 3 and 4)
    M̃ -= F[3] * (Diagonal(1 ./ (diag(τ[3]))) * HfI[3]) * F[3]'
    M̃ -= F[4] * (Diagonal(1 ./ (diag(τ[4]))) * HfI[4]) * F[4]'
  
    HfI_FT = (HfI_F1T, HfI_F2T, HfI_F3T, HfI_F4T)
    return (M̃ , F, τ, H̃, HfI_FT)
end

function bdry_vec_mod!(g, F, τ, r, s, bc_Dirichlet, bc_Neumann)

    Nr = length(r)
    Ns = length(s)

    g[:] .= 0

    # FACE 1 (Dirichlet):
    vf = bc_Dirichlet(1, 0, s)
    g[:] -= F[1] * vf

    # FACE 2 (Dirichlet):
    vf = bc_Dirichlet(2, r[Nr], s)
    g[:] -= F[2] * vf

    # FACE 3 (Neumann):
    gN = bc_Neumann(3, r, 0, 0, -1)
    vf = gN  ./ diag(τ[3])
    g[:] -= F[3] * vf

    # FACE 4 (Neumann):
    gN = bc_Neumann(4, r, s[Ns], 0, 1)
    vf = gN  ./ diag(τ[4])
    g[:] -= F[4] * vf 
    

end

function computetraction_stripped(HfI_FT, τ, lf, u, δ)
    HfI_FT = HfI_FT[lf]
    τf = τ[lf]

    return (HfI_FT * u + τf * (δ .- δ / 2)) 
  end
  

  function rateandstate(V, psi, σn, ϕ, η, a, V0)
    Y = (1 ./ (2 .* V0)) .* exp.(psi ./ a)
    f = a .* asinh.(V .* Y)
    dfdV  = a .* (1 ./ sqrt.(1 + (V .* Y).^2)) .* Y
  
    g    = σn .* f    + η .* V - ϕ
    dgdV = σn .* dfdV + η
    (g, dgdV)
  end
  
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
  