
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
    return SinModel(f, coef(fit)...)
end

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
    return SinModel(coef(fit)...)
end

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
    return MixedLinSinModel(f, coef(fit)...)
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
    return MixedLinSinModel(coef(fit)...)
end
