module SinusoidalRegressions

using LsqFit: curve_fit, coef
using RecipesBase
using PrecompileTools
using Zygote
using LinearAlgebra

# abstract types
export SRAlgorithm, SRModel, SRProblem

# main functions
export rmse, mae, torect, sinfit

# regression models
export SinModel, MixedLinSinModel

# problems
export Sin3Problem, Sin4Problem, MixedLinSin4Problem, MixedLinSin5Problem

# algorithms
export IEEE1057, IntegralEquations, LevMar, Liang

include("typedefs.jl")
include("jacquelin_sr.jl")
include("jacquelin_mlsr.jl")
include("jacquelin_gen.jl")
include("jacquelin_damped.jl")
include("ieee1057.jl")
include("lsqfit.jl")
include("liang.jl")
include("plotrecipes.jl")

"""
    sinfit(problem::Problem, algorithm::Algorithm)

Calculate a sinosoidal regression on `problem` using `algorithm`.

Currently supported problem types are:
* `Sin3Problem` -- three-parameter sinusoidal regression
* `Sin4Problem` -- four-parameter sinusoidal regression
* `MixedLinSin4Problem` -- four-parameter mixed linear and sinusoidal regression
* `MixedLinSin5Problem` -- five-parameter mixed linear and sinusoidal regression

Currently supported algorithms are:
* `IEEE1057`
* `IntegralEquations`
* `LevMar`
* `Liang`

Example
=======

```
julia> using SinusoidalRegressions
julia> t = range(0, 1, length = 100)                  # time instants
julia> s = sin.(2*pi*15*t .+ pi/4) .+ 0.1*randn(100)  # noisy samples
julia> p = Sin3Problem(t, s, 15)                      # define regression problem
julia> sinfit(p, IEEE1057())                          # calculate fit with IEEE 1057
Sinusoidal parameters SinModel{Float64}:
  Frequency (Hz)      : 15.0
  DC                  : -0.01067218324878172
  Sine amplitude (Q)  : 0.7299806464221965
  Cosine amplitude (I): 0.6822068658523716
```

See the documentation for more details.
"""
sinfit(p::SRProblem, a::SRAlgorithm) = throw("Not implemented.")

function sinfit(p::Sin3Problem, ::IEEE1057)
    ieee1057(p)
end

function sinfit(p::Sin3Problem, a::LevMar ; kwargs...)
    levmar(p, a ; kwargs...)
end

function sinfit(p::Sin4Problem, a::IEEE1057)
    ieee1057(p, a)
end

function sinfit(p::Sin4Problem, ::IntegralEquations)
    jacquelin(p)
end

function sinfit(p::Sin4Problem, a::LevMar ; kwargs...)
    levmar(p, a ; kwargs...)
end

function sinfit(p::Sin4Problem, a::Liang)
    liang(p, a)
end

function sinfit(p::MixedLinSin4Problem, a::LevMar ; kwargs...)
    levmar(p, a ; kwargs...)
end

function sinfit(p::MixedLinSin5Problem, ::IntegralEquations)
    jacquelin(p)
end

function sinfit(p::MixedLinSin5Problem, a::LevMar ; kwargs...)
    levmar(p, a ; kwargs...)
end

"""
    rmse(fit::T, exact::T, x) where {T <: SRModel}

Calculate the root mean-square error between `fit` and `exact` sampled at collection `x`.

See also: [`mae`](@ref)
"""
function rmse(fit::T, exact::T, x) where {T <: SRModel}
    fitvalues = fit(x)
    exactvalues = exact(x)
    return sqrt(sum((fitvalues .- exactvalues).^2) / length(x))
end

"""
    rmse(fit::T, samples, x) where {T <: SRModel}

Calculate the root mean-square error between `fit` and the `samples` taken at `x`.

See also: [`mae`](@ref)
"""
function rmse(fit::T, samples::AbstractVector, x) where  {T <: SRModel}
    fitvalues = fit(x)
    return sqrt(sum((fitvalues .- samples).^2) / length(x))
end

"""
    mae(fit::T, exact::T, x) where {T <: SRModel}

Calculate the mean absolute error between `fit` and `exact` sampled at collection `x`.

See also: [`rmse`](@ref)
"""
function mae(fit::T, exact::T, x) where {T <: SRModel}
    fitvalues = fit(x)
    exactvalues = exact(x)
    return sum(abs.(fitvalues .- exactvalues))/length(x)
end

"""
    torect(M, θ)

Convert polar coordinates to rectangular coordinates `(y, x)` where ``x = M\\cos(θ)`` and
``y = -M\\sin(θ)``
"""
torect(M, θ) = (-M*sin(θ), M*cos(θ))

## Precompilation
@setup_workload begin
    x = [-1.983, -1.948, -1.837, -1.827, -1.663, -0.815, -0.778, -0.754, -0.518,  0.322,  0.418,  0.781,
          0.931,  1.510,  1.607]
    y = [ 0.936,  0.810,  0.716,  0.906,  0.247, -1.513, -1.901, -1.565, -1.896,  0.051,  0.021,  1.069,
          0.862,  0.183,  0.311]
    @compile_workload begin
        p1 = Sin3Problem(x, y, 1.5)
        t2 = sinfit(p1, IEEE1057())
        t3 = sinfit(p1, LevMar())

        p2 = Sin4Problem(x, y, f = 1.5)
        t4 = sinfit(p2, IEEE1057())
        t5 = sinfit(p2, IntegralEquations())
        t6 = sinfit(p2, LevMar())
        t7 = sinfit(p2, Liang())

        p3 = MixedLinSin4Problem(x, y, 1.5)
        t8 = sinfit(p3, LevMar())

        p4 = MixedLinSin5Problem(x, y)
        t9 = sinfit(p4, IntegralEquations())
    end
end

end  # module
