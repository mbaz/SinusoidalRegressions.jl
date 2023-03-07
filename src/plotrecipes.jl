# plot recipes

"Plot a regression"
@recipe function f(t, p::T ; fitlabel = "") where {T <: SinusoidalFunctionParameters}
    xguide --> "x"
    yguide --> "y"
    ttl = ""
    if p isa SinusoidP
        ttl = "Sinusoidal Curve"
    elseif fit isa MixedLinearSinusoidP
        ttl = "Mixed Linear-Sinusoid Fit"
    end
    plot_title --> ttl
    plot_titlefontsize --> 12
    label --> fitlabel
    linecolor --> :blue
    t, p.(t)
end

# plot data, best fit, and optionally the exact function
@recipe function f(X, Y, fit::T ;
                   exact      = nothing,
                   samples    = 100,
                   fitlabel   = "fit",
                   datalabel  = "data",
                   exactlabel = "exact") where {T <: SinusoidalFunctionParameters}
    xguide --> "x"
    yguide --> "y"
    ttl = ""
    if fit isa SinusoidP
        ttl = "Sinusoidal Curve Fit"
    elseif fit isa MixedLinearSinusoidP
        ttl = "Mixed Linear-Sinusoid Fit"
    end
    plot_title --> ttl
    plot_titlefontsize --> 12
    # data (samples)
    @series begin
        seriestype --> :scatter
        label --> datalabel
        X, Y
    end
    # fit
    @series begin
        label --> fitlabel
        linecolor --> :blue
        t = range(first(X), last(X), length = samples)
        t, fit.(t)
    end
    # exact
    if !isnothing(exact)
        if exact isa T
            @series begin
                label --> exactlabel
                linecolor --> :black
                linestyle --> :dash
                t = range(first(X), last(X), length = samples)
                t, exact.(t)
            end
        else
            error("Parameters must be of the same type; got $T and $(typeof(exact))")
        end
    end
end
