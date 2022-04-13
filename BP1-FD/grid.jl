
⊗(A,B) = kron(A,B)

function transforms_e(Lw, r̂, l)
    
    A = ((Lw/2) - Lw*r̂ - Lw)/(2*tanh((r̂-1)/l) + tanh(-2/l)*(r̂ - 1))
        
    b = (A*tanh(-2/l) + Lw)/2
    c = Lw - b
        
    xt = (r,s) -> (A .* tanh.((r .- 1) ./ l) .+ b .* r .+ c,
                   ((A .* sech.((r .- 1) ./ l).^2) ./ l) .+ b,
                   zeros(size(s)))
    
    yt = (r,s) -> (A .* tanh.((s .- 1) ./ l) .+ b.*s .+ c,
                   zeros(size(r)),
                   ((A .* sech.((s .- 1) ./ l).^2) ./ l) .+ b)
    
    return xt, yt
    
end


function create_metrics(ne, μ,
                        xf=(r,s)->(r, ones(size(r)), zeros(size(r))),
                        yf=(r,s)->(s, zeros(size(s)), ones(size(s))))
    nn = ne + 1
    N = nn^2
    

    r = range(-1, stop=1, length=nn)
    s = range(-1, stop=1, length=nn)

    # Create the mesh
    r = ones(1, nn) ⊗ r
    s = s' ⊗ ones(nn)
    (x, xr, xs) = xf(r, s)
    (y, yr, ys) = yf(r, s)

    J = xr .* ys - xs .* yr
    @assert minimum(J) > 0
    JI = 1 ./ J
    
    rx =  ys ./ J
    sx = -yr ./ J
    ry = -xs ./ J
    sy =  xr ./ J
    
    crr = J .* (rx .* μ .* rx + ry .* μ .* ry)
    crs = J .* (sx .* μ .* rx + sy .* μ .* ry)
    css = J .* (sx .* μ .* sx + sy .* μ .* sy)
    
    #
    # Block surface matrices
    #
    (xf1, yf1) = (view(x, 1, :), view(y, 1, :))
    nx1 = -ys[1, :]
    ny1 =  xs[1, :]
    sJ1 = hypot.(nx1, ny1)
    nx1 = nx1 ./ sJ1
    ny1 = ny1 ./ sJ1

    (xf2, yf2) = (view(x, nn, :), view(y, nn, :))
    nx2 =  ys[end, :]
    ny2 = -xs[end, :]
    sJ2 = hypot.(nx2, ny2)
    nx2 = nx2 ./ sJ2
    ny2 = ny2 ./ sJ2

    (xf3, yf3) = (view(x, :, 1), view(y, :, 1))
    nx3 =  yr[:, 1]
    ny3 = -xr[:, 1]
    sJ3 = hypot.(nx3, ny3)
    nx3 = nx3 ./ sJ3
    ny3 = ny3 ./ sJ3

    (xf4, yf4) = (view(x, :, nn), view(y, :, nn))
    nx4 = -yr[:, end]
    ny4 =  xr[:, end]
    sJ4 = hypot.(nx4, ny4)
    nx4 = nx4 ./ sJ4
    ny4 = ny4 ./ sJ4


    (coord = (x,y),
     facecoord = ((xf1, xf2, xf3, xf4), (yf1, yf2, yf3, yf4)),
     crr = crr, css = css, crs = crs,
     J=J,
     JI = JI,
     sJ = (sJ1, sJ2, sJ3, sJ4),
     nx = (nx1, nx2, nx3, nx4),
     ny = (ny1, ny2, ny3, ny4),
     rx = rx, ry = ry, sx = sx, sy = sy)
end
