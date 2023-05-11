
function jacquelin(prob::MixedLinSin5Problem)
    (; X, Y) = prob
    # verify elements of X are ordered
    if !issorted(X)
        error("Please provide abscissa vector in ascending order.")
    end

    _jacquelin_mixlinsin(X, Y)
end

function _jacquelin_mixlinsin(X, Y)
    n = length(X)

    # "short way" -- actually, a first approximation
    S = zeros(eltype(X), n)
    for k = 2:n
        S[k] = S[k-1] + 0.5*(Y[k]+Y[k-1])*(X[k]-X[k-1])
    end

    SS = zeros(eltype(X), n)
    for k = 2:n
        SS[k] = SS[k-1] + 0.5*(S[k]+S[k-1])*(X[k]-X[k-1])
    end

    ΣSS² = sum(SS.^2)
    Σx³SS = sum(SS .* X.^3)
    Σx²SS = sum(SS .* X.^2)
    ΣxSS = sum(SS .* X)
    ΣSS = sum(SS)
    Σx = sum(X)
    Σx² = sum(X.^2)
    Σx³ = sum(X.^3)
    Σx⁴ = sum(X.^4)
    Σx⁵ = sum(X.^5)
    Σx⁶ = sum(X.^6)
    Σy = sum(Y)
    Σyx = sum(X.*Y)
    Σyx² = sum(X.^2 .* Y)
    Σyx³ = sum(X.^3 .* Y)
    ΣySS = sum(Y .* SS)

    M = [ ΣSS²  Σx³SS Σx²SS ΣxSS ΣSS ;
          Σx³SS Σx⁶   Σx⁵   Σx⁴  Σx³ ;
          Σx²SS Σx⁵   Σx⁴   Σx³  Σx² ;
          ΣxSS  Σx⁴   Σx³   Σx²  Σx  ;
          ΣSS   Σx³   Σx²   Σx   n   ]

    x = [ΣySS, Σyx³, Σyx², Σyx, Σy]

    (A0, E0, B0, C0, D0) = M \ x

    ω1 = sqrt(-A0)
    β = sin.(ω1.*X)
    η = cos.(ω1.*X)

    Σβ = sum(β)
    Ση = sum(η)
    Σxβ = sum(X .* β)
    Σxη = sum(X .* η)
    Σβ² = sum(β.^2)
    Ση² = sum(η.^2)
    Σβη = sum(β .* η)
    Σyβ = sum(β .* Y)
    Σyη = sum(η .* Y)

    N = [ n  Σx  Σβ  Ση ;
          Σx Σx² Σxβ Σxη ;
          Σβ Σxβ Σβ² Σβη ;
          Ση Σxη Σβη Ση² ]

    x = [Σy, Σyx, Σyβ, Σyη]

    (a1, p1, b1, c1) = N \ x

    ρ1 = sqrt(b1^2 + c1^2)
    φ1 = b1 > 0 ? atan(c1/b1) : atan(c1/b1) + π

    # "full way" -- finer approximation
    K = round.( (ω1*X .+ φ1)/π )
    θ = zeros(eltype(X), n)

    for k = 1:n
        r = Y[k] - a1 - p1*X[k]
        if ρ1^2 > r^2
            θ[k] = (-1)^K[k] * atan(r / sqrt(ρ1^2 - r^2)) + K[k]*π
        else
            if r > 0
                θ[k] = (-1)^K[k] * π/2 + K[k]*π
            else
                θ[k] = -(-1)^K[k] * π/2 + K[k]*π
            end
        end
    end

    O = [ Σx² Σx ;
          Σx  n  ]

    x = [sum(θ .* X), sum(θ)]

    (ω2, φ2) = O \ x

    β = sin.(ω2.*X)
    η = cos.(ω2.*X)

    Σβ = sum(β)
    Ση = sum(η)
    Σxβ = sum(X .* β)
    Σxη = sum(X .* η)
    Σβ² = sum(β.^2)
    Ση² = sum(η.^2)
    Σβη = sum(β .* η)
    Σyβ = sum(β .* Y)
    Σyη = sum(η .* Y)

    P = [ n  Σx  Σβ  Ση ;
          Σx Σx² Σxβ Σxη ;
          Σβ Σxβ Σβ² Σβη ;
          Ση Σxη Σβη Ση² ]

    x = [Σy, Σyx, Σyβ, Σyη]

    (a3, p3, b3, c3) = P \ x

    return MixedLinSinModel(ω2/(2π), a3, b3, c3, p3)

end
