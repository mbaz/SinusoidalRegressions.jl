"""
    sinfit(X, Y, guess::SinusoidP = sinfit_j(X, Y) ; kwargs...)

Perform a four-parameter least-squares sinusoidal fit of the independent variables
`X` and dependent variables `Y`, using the non-linear optimization solver from
`LsqFit.jl`.

The data is fit to the model ``s(x; f, DC, Q, I) = DC + Q\\sin(2πfx) + I\\cos(2πfx)``. The
values in `X` must be sorted in ascending order.

The Levenberg-Marquardt algorithm used by `LsqFit.jl` requires an initial guess
of the parameters `f`, `DC`, `Q` and `I`. If no initial guess is provided, then
one is calculated using [`sinfit_j`](@ref).

All keyword arguments provided are directly passed to `LsqFit.curve_fit`.

See also [`SinusoidP`](@ref), [`sinfit_j`](@ref), [LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl),
[`curve_fit`](https://github.com/JuliaNLSolvers/LsqFit.jl)
"""
function sinfit(X, Y, guess::SinusoidP = sinfit_j(X, Y) ; kwargs...)
    # Convert guess from a SP4 struct to an array
    guessarr = [guess.f, guess.DC, guess.Q, guess.I]
    # Define data model
    @. model(t, p) = p[2] + p[3]*sin(2π*p[1]*t) + p[4]*cos(2π*p[1]*t)
    # fit
    fit = curve_fit(model, X, Y, guessarr ; kwargs...)
    return SinusoidP(coef(fit)...)
end

"""
    sinfit(X, Y, f ; kwargs...)

Perform a three-parameter least-squares sinusoidal fit of the independent variables
`X` and dependent variables `Y`, assuming a known frequency `f`, using the non-linear
optimization solver from `LsqFit.jl`.

The data is fit to the model ``s(x; DC, Q, I) = DC + Q\\sin(2πfx) + I\\cos(2πfx)``. The
values in `X` must be sorted in ascending order.

The Levenberg-Marquardt algorithm used by `LsqFit.jl` requires an initial guess
of the parameters, `DC`, `Q` and `I`. If no initial guess is provided, then one is
calculated using the linear regression algorithm from IEEE 1057 (see [`ieee1057`](@ref)).

All keyword arguments provided are directly passed to `LsqFit.curve_fit`.

See also [`SinusoidP`](@ref), [`ieee1057`](@ref), [`curve_fit`](https://github.com/JuliaNLSolvers/LsqFit.jl).
"""
function sinfit(X, Y, f ; kwargs...)
    # pre-fit with IEEE 1057
    (; DC, Q, I) = _ieee1057(X, Y, f)
    # define data model
    guessarr = [DC, Q, I]
    @. model(t, p) = p[1] + p[2]*sin(2π*f*t) + p[3]*cos(2π*f*t)
    # fit
    fit = curve_fit(model, X, Y, guessarr ; kwargs...)
    return SinusoidP(f, coef(fit)...)
end

"""
    mixlinsinfit(X, Y, guess::MixedLinearSinusoidP = mixlinsinfit_j(X, Y) ; kwargs...)

Perform a least-squares mixed linear-sinusoid fit of the independent variables `X` and dependent
variables `Y`,  using the non-linear optimization solver from `LsqFit.jl`.

The data is fit to the model ``s(x, f, DC, Q, I, m) = DC + Q\\sin(2πfx) + I\\cos(2πfx) + mx``. The
values in `X` must be sorted in axcending order.

The Levenberg-Marquardt algorithm used by `LsqFit.jl` requires an initial guess
of the parameters, `DC`, `Q` `I`, and `m`. If no initial guess is provided, then one is
calculated then one is calculated using [`mixlinsinfit_j`](@ref).

All keyword arguments provided are directly passed to `LsqFit.curve_fit`.

See also [`MixedLinearSinusoidP`](@ref), [`mixlinsinfit_j`](@ref),
[`curve_fit`](https://github.com/JuliaNLSolvers/LsqFit.jl).
"""
function mixlinsinfit(X, Y, guess::MixedLinearSinusoidP = mixlinsinfit_j(X, Y) ; kwargs...)
    # Convert guess to array
    guessarr = [guess.f, guess.DC, guess.Q, guess.I, guess.m]
    # Define data model
    @. model(t, p) = p[2] + p[3]*sin(2π*p[1]*t) + p[4]*cos(2π*p[1]*t) + t*p[5]
    # fit
    fit = curve_fit(model, X, Y, guessarr ; kwargs...)
    return MixedLinearSinusoidP(coef(fit)...)
end
