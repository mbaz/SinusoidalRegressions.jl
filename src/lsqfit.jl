"""
    sinfit(X, Y, f ; kwargs...)

Perform a three-parameter least-squares sinusoidal fit of the independent variables
`X` and dependent variables `Y`, assuming a known frequency `f`, using the non-linear
optimization solver from `LsqFit.jl`.

The data is fit to the model ``s(x; DC, Q, I) = DC + Q\\sin(2πfx) + I\\cos(2πfx)``. The
values in `X` must be sorted in ascending order.

The Levenberg-Marquardt algorithm used by `LsqFit.jl` requires an initial guess
of all three parameters, `DC`, `Q` and `I`. If no initial guess is provided, then one is
calculated using the linear regression algorithm from IEEE 1057 (see [`ieee1057`](@ref)).

All keyword arguments provided are directly passed to `LsqFit.curve_fit`.

See also [`SinusoidP`](@ref), [`ieee1057`](@ref),
[`curve_fit`](https://github.com/JuliaNLSolvers/LsqFit.jl).
"""
function levmar(prob::Sin3Problem, a::LevMar ; kwargs...)
    (; X, Y, f, DC, Q, I, lb, ub) = prob
    # Check if all parameters are given initial estimates
    if any(map(ismissing, (DC, Q, I)))
        # if not all guesses are given, find our own estimates
        initialfit = _jacquelin_sin(X, Y)
        ismissing(DC) && (DC = initialfit.DC)
        ismissing(Q)  && (Q  = initialfit.Q)
        ismissing(I)  && (I  = initialfit.I)
    end
    guess = [DC, Q, I]
    # bounds
    ismissing(lb) && (lb = [-Inf, -Inf, -Inf])
    ismissing(ub) && (ub = [+Inf, +Inf, +Inf])
    # Define data model
    @. model(t, p) = p[1] + p[2]*sin(2π*f*t) + p[3]*cos(2π*f*t)
    # fit
    if a.use_ga
        function Avv!(dir_deriv,p,v)
            for i = 1:length(X)
                dir_deriv[i] = transpose(v) * Zygote.hessian(z -> model(X[i], z), p) * v
            end
        end
        # fit
        fit = curve_fit(model, X, Y, guess,
                        avv! = Avv!, lambda = 0, min_step_quality = 0,
                        lower = lb, upper = ub, kwargs...)
    else
        fit = curve_fit(model, X, Y, guess, lower = lb, upper = ub, kwargs...)
    end
    return SinusoidP(f, coef(fit)...)
end

"""
    sinfit(X, Y, guess::SinusoidP = sinfit_j(X, Y) ; kwargs...)

Perform a four-parameter least-squares sinusoidal fit of the independent variables
`X` and dependent variables `Y`, using the non-linear optimization solver from
`LsqFit.jl`.

The data is fit to the model ``s(x; f, DC, Q, I) = DC + Q\\sin(2πfx) + I\\cos(2πfx)``. The
values in `X` must be sorted in ascending order.

The Levenberg-Marquardt algorithm used by `LsqFit.jl` requires an initial guess
of all four model parameters. If no initial guess is provided, then one is calculated
using [`sinfit_j`](@ref).

All keyword arguments provided are directly passed to `LsqFit.curve_fit`.

See also [`SinusoidP`](@ref), [`sinfit_j`](@ref),
[LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl),
[`curve_fit`](https://github.com/JuliaNLSolvers/LsqFit.jl)
"""
function levmar(prob::Sin4Problem, a::LevMar ; kwargs...)
    (; X, Y, f, DC, Q, I, lb, ub) = prob
    # Check if all parameters are given initial estimates
    if any(map(ismissing, (f, DC, Q, I)))
        # if not all guesses are given, find our own estimates
        initialfit = _jacquelin_sin(X, Y)
        ismissing(f)  && (f  = initialfit.f)
        ismissing(DC) && (DC = initialfit.DC)
        ismissing(Q)  && (Q  = initialfit.Q)
        ismissing(I)  && (I  = initialfit.I)
    end
    guess = [f, DC, Q, I]
    # bounds
    ismissing(lb) && (lb = [-Inf, -Inf, -Inf, -Inf])
    ismissing(ub) && (ub = [+Inf, +Inf, +Inf, +Inf])
    # define data model
    guessarr = [f, DC, Q, I]
    @. model(t, p) = p[2] + p[3]*sin(2π*p[1]*t) + p[4]*cos(2π*p[1]*t)
    # fit
    if a.use_ga
        function Avv!(dir_deriv,p,v)
            for i = 1:length(X)
                dir_deriv[i] = transpose(v) * Zygote.hessian(z -> model(X[i], z), p) * v
            end
        end
        fit = curve_fit(model, X, Y, guess,
                        avv! = Avv!, lambda = 0, min_step_quality = 0,
                        lower = lb, upper = ub, kwargs...)
    else
        fit = curve_fit(model, X, Y, guessarr ; kwargs...)
    end
    return SinusoidP(coef(fit)...)
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
function levmar(prob::MixedLinSin4Problem, a::LevMar ; kwargs...)
    (; X, Y, f, DC, Q, I, m, lb, ub) = prob
    # Check if all parameters are given initial estimates
    if any(map(ismissing, (DC, Q, I, m)))
        # if not all guesses are given, find our own estimates
        initialfit = _jacquelin_mixlinsin(X, Y)
        ismissing(DC) && (DC = initialfit.DC)
        ismissing(Q)  && (Q  = initialfit.Q)
        ismissing(I)  && (I  = initialfit.I)
        ismissing(m)  && (m  = initialfit.m)
    end
    guess = [DC, Q, I, m]
    # bounds
    ismissing(lb) && (lb = [-Inf, -Inf, -Inf, -Inf])
    ismissing(ub) && (ub = [+Inf, +Inf, +Inf, +Inf])
    # Define data model
    @. model(t, p) = p[1] + p[2]*sin(2π*f*t) + p[3]*cos(2π*f*t) + t*p[4]
    # fit
    if a.use_ga
        function Avv!(dir_deriv,p,v)
            for i = 1:length(X)
                dir_deriv[i] = transpose(v) * Zygote.hessian(z -> model(X[i], z), p) * v
            end
        end
        # fit
        fit = curve_fit(model, X, Y, guess,
                        avv! = Avv!, lambda = 0, min_step_quality = 0,
                        lower = lb, upper = ub, kwargs...)
    else
        fit = curve_fit(model, X, Y, guess, lower = lb, upper = ub, kwargs...)
    end
    return MixedLinearSinusoidP(f, coef(fit)...)
end

function levmar(prob::MixedLinSin5Problem, a::LevMar ; kwargs...)
    (; X, Y, f, DC, Q, I, m, lb, ub) = prob
    # Check if all parameters are given initial estimates
    if any(map(ismissing, (f, DC, Q, I, m)))
        # if not all guesses are given, find our own estimates
        initialfit = _jacquelin_mixlinsin(X, Y)
        ismissing(f)  && (f  = initialfit.f)
        ismissing(DC) && (DC = initialfit.DC)
        ismissing(Q)  && (Q  = initialfit.Q)
        ismissing(I)  && (I  = initialfit.I)
        ismissing(m)  && (m  = initialfit.m)
    end
    guess = [f, DC, Q, I, m]
    # bounds
    ismissing(lb) && (lb = [0, -Inf, -Inf, -Inf, -Inf])
    ismissing(ub) && (ub = [+Inf, +Inf, +Inf, +Inf, +Inf])
    # Define data model
    @. model(t, p) = p[2] + p[3]*sin(2π*p[1]*t) + p[4]*cos(2π*p[1]*t) + t*p[5]
    # fit
    if a.use_ga
        function Avv!(dir_deriv,p,v)
            for i = 1:length(X)
                dir_deriv[i] = transpose(v) * Zygote.hessian(z -> model(X[i], z), p) * v
            end
        end
        # fit
        fit = curve_fit(model, X, Y, guess,
                        avv! = Avv!, lambda = 0, min_step_quality = 0,
                        lower = lb, upper = ub, kwargs...)
    else
        fit = curve_fit(model, X, Y, guess, lower = lb, upper = ub, kwargs...)
    end
    return MixedLinearSinusoidP(coef(fit)...)
end
