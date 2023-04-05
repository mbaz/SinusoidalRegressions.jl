module SinusoidalRegressions

using LsqFit
using RecipesBase
using SnoopPrecompile

export sinfit_j, mixlinsinfit_j,
       sinfit, mixlinsinfit,
       ieee1057,
       rmse, mae,
       SinusoidalFunctionParameters,
       SinusoidP, MixedLinearSinusoidP,
       torect

include("typedefs.jl")
include("jacquelin_sr.jl")
include("jacquelin_mlsr.jl")
include("jacquelin_gen.jl")
include("jacquelin_damped.jl")
include("ieee1057.jl")
include("lsqfit.jl")
include("plotrecipes.jl")

"""
    rmse(fit::T, exact::T, x) where {T <: SinusoidalFunctionParameters}

Calculate the root mean-square error between `fit` and `exact` sampled at collection `x`.

See also: [`mae`](@ref)
"""
function rmse(fit::T, exact::T, x) where {T <: SinusoidalFunctionParameters}
    fitvalues = fit(x)
    exactvalues = exact(x)
    return sqrt(sum((fitvalues .- exactvalues).^2) / length(x))
end

"""
    mae(fit::T, exact::T, x) where {T <: SinusoidalFunctionParameters}

Calculate the mean absolute error between `fit` and `exact` sampled at collection `x`.

See also: [`rmse`](@ref)
"""
function mae(fit::T, exact::T, x) where {T <: SinusoidalFunctionParameters}
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
@precompile_setup begin
    x = [-1.983, -1.948, -1.837, -1.827, -1.663, -0.815, -0.778, -0.754, -0.518,  0.322,  0.418,  0.781,
          0.931,  1.510,  1.607]
    y = [ 0.936,  0.810,  0.716,  0.906,  0.247, -1.513, -1.901, -1.565, -1.896,  0.051,  0.021,  1.069,
          0.862,  0.183,  0.311] 
    @precompile_all_calls begin
        t2 = ieee1057(x, y)
        t1 = ieee1057(x, y, 1.5)
        t3 = sinfit_j(x, y)
        t4 = sinfit(x, y, t3)
        t5 = mixlinsinfit_j(x, y)
        t6 = mixlinsinfit(x, y, t5)
    end
end

end  # module
