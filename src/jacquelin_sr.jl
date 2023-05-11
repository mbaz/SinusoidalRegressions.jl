
function jacquelin(p::Sin4Problem)
    (; X, Y) = p
    n = length(X)
    return _jacquelin_sin(X, Y)
end

function _jacquelin_sin(X, Y)
    # verify elements of X are ordered
    if !issorted(X)
        error("Please provide abscissa vector in ascending order.")
    end

    sec1 = _jacquelin_part1(X, Y)       # first section
    sec2 = _jacquelin_part2(X, Y, sec1) # second section
    sec3 = _jacquelin_part3(X, Y, sec2) # third section
end

function _jacquelin_part1(X, Y)
    n = length(X)

    S = zeros(eltype(X), n)
    for k = 2:n
        S[k] = S[k-1] + 0.5*(Y[k]+Y[k-1])*(X[k]-X[k-1])
    end

    SS = zeros(eltype(X), n)
    for k = 2:n
        SS[k] = SS[k-1] + 0.5*(S[k]+S[k-1])*(X[k]-X[k-1])
    end

    Σx = sum(X)
    Σx² = sum(X.^2)
    Σx³ = sum(X.^3)
    Σx⁴ = sum(X.^4)
    ΣSS = sum(SS)
    ΣSS² = sum(SS.^2)
    ΣxSS = sum(X .* SS)
    Σx²SS = sum(X.^2 .* SS)
    Σy = sum(Y)
    Σyx = sum(X .* Y)
    Σyx² = sum(X.^2 .* Y)
    ΣySS = sum(Y .* SS)

    M = [ ΣSS²   Σx²SS  ΣxSS  ΣSS
          Σx²SS  Σx⁴    Σx³   Σx²
          ΣxSS   Σx³    Σx²   Σx
          ΣSS    Σx²    Σx    n   ]

    (A1, B1, C1, D1) = M \ [ΣySS, Σyx², Σyx, Σy]

    x1 = X[1]
    ω1 = sqrt(-A1)
    a1 = 2*B1 / (ω1^2)
    b1 = (B1*x1^2 + C1*x1 + D1 - a1) * sin(ω1*x1) + (1/ω1)*(C1 + 2*B1*x1) * cos(ω1*x1)
    c1 = (B1*x1^2 + C1*x1 + D1 - a1) * cos(ω1*x1) - (1/ω1)*(C1 + 2*B1*x1) * sin(ω1*x1)

    return SinModel(ω1/(2π), a1, b1, c1) 
end

function _jacquelin_part2(X, Y, sp)
    n = length(X)
    ω1, a1, b1, c1 = 2π*sp.f, sp.DC, sp.Q, sp.I

    a2 = a1
    ρ2 = sqrt(b1^2 + c1^2)
    if b1 == 0
        φ1 = π/2 * sign(c1)
    else
        if b1 > 0
            φ1 = atan(c1/b1)
        else
            φ1 = atan(c1/b1) + π
        end
    end

    K = [round((ω1*xk + φ1)/π) for xk in X]

    θ = zeros(n)
    for k = 1:n
        if ρ2^2 > (Y[k] - a2)^2
            θ[k] = (-1)^K[k] * atan( (Y[k]-a2) / sqrt(ρ2^2 - (Y[k]-a2)^2) ) + π*K[k]
        else
            if Y[k] > a2
                θ[k] = (-1)^K[k]*π/2 + K[k]*π
            else
                θ[k] = -(-1)^K[k]*π/2 + K[k]*π
            end
        end
    end

    Σx  = sum(X)
    Σx² = sum(X.^2)
    Σθ  = sum(θ)
    Σθx = sum(θ .* X)

    M = [ Σx²  Σx
          Σx   n  ]

    (ω2, φ2) = M \ [ Σθx, Σθ ]

    b2 = ρ2 * cos(φ2)
    c2 = ρ2 * sin(φ2)

    return SinModel(ω2/(2π), a2, b2, c2)
end

function _jacquelin_part3(X, Y, sp)
    n = length(X)
    ω3 = 2π*sp.f

    Σsin = sum(sin, ω3*X)
    Σcos = sum(cos, ω3*X)
    Σsin² = sum(x -> sin(x)*sin(x), ω3*X)
    Σcos² = sum(x -> cos(x)*cos(x), ω3*X)
    Σsincos = sum(x -> sin(x)*cos(x), ω3*X)
    Σy = sum(Y)
    Σysin = sum( (Y[k]*sin(ω3*X[k]) for k in eachindex(X)) )
    Σycos = sum( (Y[k]*cos(ω3*X[k]) for k in eachindex(X)) )

    M = [ n     Σsin     Σcos
          Σsin  Σsin²    Σsincos
          Σcos  Σsincos  Σcos²   ]

    (a3, b3, c3) = M \ [ Σy, Σysin, Σycos ]

    return SinModel(ω3/(2π), a3, b3, c3)
end
