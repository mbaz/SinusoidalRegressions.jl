
function ieee1057(p::Sin3Problem)
    (; X, Y, f) = p
    if length(X) != length(Y)
        error("Input vectors must be of equal length (got $(length(X)) and $(length(Y))).")
    end
    if f <= 0
        error("Frequency must be larger than 0.")
    end
    return _ieee1057_3P(X, Y, f)
end

function ieee1057(p::Sin4Problem, a::IEEE1057)
    (; X, Y, f) = p
    if length(X) != length(Y)
        error("Input vectors must be of equal length.")
    end
    if !ismissing(f) && f < 0
        error("Frequency estimate must be larger than 0.")
    end
    return _ieee1057_4P(X, Y, f, a.iterations)
end

# 3-parameter fit
function _ieee1057_3P(X, Y, f)
    M = length(X)

    D0 = hcat([1 for k in eachindex(X)],
              [sin(2*π*f*x) for x in X],
              [cos(2*π*f*x) for x in X])

    (DC, Q, I) = D0 \ Y

    return SinModel(f, DC, Q, I)
end

# 4-parameter fit
function _ieee1057_4P(X, Y, f, iterations)
    if ismissing(f)
        # obtain a frequency estimate using Jacquelin's regression
        sec1 = _jacquelin_part1(X, Y)
        f = _jacquelin_part2(X, Y, sec1).f
    end

    # pre-fit using 3-parameter IEEE 1057 algorithm
    (;DC, Q, I) = _ieee1057_3P(X, Y, f)

    # Pre-allocations
    D = zeros(eltype(X), length(X), 4)
    D[:, 3] .= 1.0
    Δf = 0.0

    # iteration
    for i = 1:iterations
        f += Δf
        D[:, 1] .= cos.(2π*f*X)
        D[:, 2] .= sin.(2π*f*X)
        D[:, 4] .= X.*(-I*sin.(2π*f*X) .+ Q*cos.(2π*f*X))

        (I, Q, DC, Δf) = D \ Y
    end
    f += Δf

    return SinModel(f, DC, Q, I)
end

"""
    ieee1057_testconvergence(X, Y ; f = nothing, iterations = 6)

Returns `true` if the 4-paramter IEEE 1057 algorithm converges, and `false` otherwise.
Convergence is determined rather crudely, by observing whether the frequency correction
Δf becomes smaller with each iteration. Argument `f` is an estimate of the sinusoidal's
frequency; if not provided, it is calculated as in [`ieee1057`](@ref).
"""
function ieee1057_testconvergence(p::Sin4Problem, a::IEEE1057 = IEEE1057())
    (; X, Y, f) = p
    if ismissing(f)
        sec1 = _jacquelin_part1(X, Y)
        f = _jacquelin_part2(X, Y, sec1).f
    end
    (;DC, Q, I) = _ieee1057_3P(X, Y, f)

    D = zeros(eltype(X), length(X), 4)
    D[:, 3] .= 1.0
    Δf = 0.0

    for i = 1:a.iterations
        prev_Δf = Δf
        f += Δf
        D[:, 1] .= cos.(2π*f*X)
        D[:, 2] .= sin.(2π*f*X)
        D[:, 4] .= X.*(-I*sin.(2π*f*X) .+ Q*cos.(2π*f*X))
        (I, Q, DC, Δf) = D \ Y
        @debug "Iteration: $i" prev_Δf Δf
        if (i > 3) && (abs(Δf) > 1e-3) && (abs(Δf) > abs(prev_Δf))
            return false
        end
    end
    return true
end
