
#
# Algorithms
#

"""
    SRAlgorithm

Supported algorithm types are subtypes of the abstract type `SRAlgorithm`.
"""
abstract type SRAlgorithm end

"""
    IEEE1057([iterations = 6]) <: SRAlgorithm

Define an instance of the IEEE 1057 sinusoidal fitting algorithm.

Optional argument `iterations` specifies how many iterations to run before the algorithm
stops. The default value is 6, which is the value recommended by the standard. This value
is only used when calculating a 4-parameter fit.
"""
Base.@kwdef struct IEEE1057 <: SRAlgorithm
    iterations::Int = 6
end

"""
    IntegralEquations() <: SRAlgorithm

Define an instance of the integral-equations sinusoidal fitting algorithm described by
J.  Jacquelin in "Régressions et équations intégrales", 2014 (available at
https://fr.scribd.com/doc/14674814/Regressions-et-equations-integrales).

This algorithm does not accept any configuration parameters.
"""
struct IntegralEquations <: SRAlgorithm end

"""
    LevMar([use_ga = false]) <: SRAlgorithm

Define an instance of the Levenberg-Marquardt sinusoidal fitting algorithm.

If the optional argument `use_ga` is set to `true`, the algorithm will use geodesic
acceleration to potentially improve its performance and accuracy. See
[`curve_fit`](https://github.com/JuliaNLSolvers/LsqFit.jl) for more details.
"""
Base.@kwdef struct LevMar <: SRAlgorithm
    use_ga :: Bool = false
end

"""
    Liang([threshold = 0.15, iterations = 100, q = 1e-5]) <: SRAlgorithm

Define an instance of the sinusoidal fitting algorithm described in Liang et al,
"Fitting Algorithm of Sine Wave with Partial Period Waveforms and Non-Uniform Sampling
Based on Least-Square Method." Journal of Physics: Conference Series 1149.1
(2018)ProQuest. Web. 17 Apr. 2023.

This algorithm is designed for scenarios where only a fraction of a period of the sinusoid
has been sampled. Its optional parameters `threshold` and `q` are described in the paper.
Additionally, a maximum number of iterations may be specified in `iterations`.
"""
Base.@kwdef struct Liang <: SRAlgorithm
    threshold  :: Float64 = 0.15
    iterations :: Int     = 100
    q          :: Float64 = 1e-5
end

#
# Problems
#

"""
    SRProblem

Supported problem types are subtypes of the abstract type `SRProblem`.
"""
abstract type SRProblem end

const MR = Union{Missing, Real}

"""
    Sin3Problem(X, Y, f , [DC, Q, I, lb, ub]) <: SRProblem

Define a three-parameter sinusoidal regression problem.

The data is fit to the model ``f(x; DC, Q, I) = DC + Q\\sin(2πfx) + I\\cos(2πfx)``. The sampling
instants are given by `X`, and the samples by `Y`. The frequency `f` (in Hz) is assumed to
be known exactly.

When using a fitting algorithm that accepts initial parameter estimates, these may be given
by the optional keyword arguments `DC`, `Q` and `I` (default `missing`). Lower and upper bounds
may be specified in `lb` and `ub`, which must be vectors of length 3.

See also: [`Sin4Problem`](@ref)
"""
Base.@kwdef struct Sin3Problem{T1, T2, F<:Real, P1<:MR, P2<:MR, P3<:MR, LB, UB} <: SRProblem
    X  :: T1  # sampling times
    Y  :: T2  # sample values
    f  :: F  # exact known frequency
    # initial estimates (may be missing)
    DC :: P1 = missing
    Q  :: P2 = missing
    I  :: P3 = missing
    lb :: LB = missing  # lower bounds
    ub :: UB = missing  # upper bounds
end

Sin3Problem(X, Y, f ; kwargs...) = Sin3Problem(; X, Y, f, kwargs...)

"""
    Sin4Problem(X, Y ; [f, DC, Q, I, lb, ub]) <: SRProblem

Define a four-parameter sinusoidal regression problem.

The data is fit to the model ``f(x; f, DC, Q, I) = DC + Q\\sin(2πfx) + I\\cos(2πfx)``.
The sampling instants are given by `X`, and the samples by `Y`.

When using a fitting algorithm that accepts initial parameter estimates, these may be given
by the optional keyword arguments `f`, `DC`, `Q` and `I` (default `missing`). Lower and upper
bounds may be specified in `lb` and `ub`, which must be vectors of length 4.

See also: [`Sin3Problem`](@ref)
"""
Base.@kwdef struct Sin4Problem{T1, T2, P1<:MR, P2<:MR, P3<:MR, P4<:MR, LB, UB} <: SRProblem
    X  :: T1  # sampling times
    Y  :: T2  # sample values
    f  :: P1 = missing  # initial estimates (may be missing)
    DC :: P2 = missing
    Q  :: P3 = missing
    I  :: P4 = missing
    lb :: LB = missing  # lower bounds
    ub :: UB = missing  # uppper bounds
end

Sin4Problem(X, Y ; kwargs...) = Sin4Problem(; X, Y, kwargs...)

"""
    MixedLinSin4Problem(X, Y, f ; [DC, Q, I, m, lb, ub]) <: SRProblem

Define a four-parameter mixed linear-sinusoidal regression problem.

The data is fit to the model ``f(x; DC, Q, I, m) = DC + Q\\sin(2πfx) + I\\cos(2πfx) + mx``.
The sampling instants are given by `X`, and the samples by `Y`. The frequency `f` (in Hz) is
assumed to be known exactly.

When using a fitting algorithm that accepts initial parameter estimates, these may be given
by the optional keyword arguments `DC`, `Q` `I` and `m` (default `missing`). Lower and upper
bounds may be specified in `lb` and `ub`, which must be vectors of length 4.

See also: [`MixedLinSin5Problem`](@ref)
"""
Base.@kwdef struct MixedLinSin4Problem{T1, T2, F<:Real, P1<:MR, P2<:MR, P3<:MR, P4<:MR, LB, UB} <: SRProblem
    X :: T1
    Y :: T2
    f :: F
    DC :: P1= missing
    Q  :: P2 = missing
    I  :: P3 = missing
    m  :: P4 = missing
    lb :: LB = missing
    ub :: UB = missing
end

MixedLinSin4Problem(X, Y, f ; kwargs...) = MixedLinSin4Problem(; X, Y, f, kwargs...)

"""
    MixedLinSin5Problem(X, Y ; [f, DC, Q, I, m, lb, ub]) <: SRProblem

Define a five-parameter mixed linear-sinusoidal regression problem.

The data is fit to the model ``f(x ; f, DC, Q, I, m) = DC + Q\\sin(2πfx) + I\\cos(2πfx) + mx``.
The sampling instants are given by `X`, and the samples by `Y`.

When using a fitting algorithm that accepts initial parameter estimates, these may be given
by the optional keyword arguments `f`, `DC`, `Q` `I` and `m` (default `missing`). Lower and upper
bounds may be specified in `lb` and `ub`, which must be vectors of length 5.

See also: [`MixedLinSin4Problem`](@ref)
"""
Base.@kwdef struct MixedLinSin5Problem{T1, T2, P1<:MR, P2<:MR, P3<:MR, P4<:MR, P5<:MR, LB, UB} <: SRProblem
    X :: T1
    Y :: T2
    f  :: P1  = missing
    DC :: P2 = missing
    Q  :: P3 = missing
    I  :: P4 = missing
    m  :: P5 = missing
    lb :: LB = missing
    ub :: UB = missing
end

MixedLinSin5Problem(X, Y ; kwargs...) = MixedLinSin5Problem(; X, Y, kwargs...)

"""
    SRModel

Supported models types are subtypes of the abstract type `SRModel`.

Supertype for the parameters of the different kinds of sinusoidal models
supported by SinusoidalRegressions.

See also: [`SinModel`](@ref), [`MixedLinSinModel`](@ref)
"""
abstract type SRModel end

"""
    SinModel{T <: Real} <: SRModel

Parameters `f`, `DC`, `Q` and `I` of the sinusoidal model ``s(x) = DC + I\\cos(2πfx) + Q\\sin(2πfx)``.

See also: [`MixedLinSinModel`](@ref).
"""
struct SinModel{T <: Real} <: SRModel
    f  :: T  # frequency in Hz
    DC :: T
    Q  :: T  # sine amplitude
    I  :: T  # cosine amplitude
end

"""
    SinModel{T}(f, DC, Q, I)

Construct a sinusoidal model with the given parameters.

To express the model as ``s(x) = M\\cos(2πfx + θ)``, use the `torect` function.

Examples
========

```
julia> SinModel(10, 1, -0.5, 1.2)
Sinusoidal parameters SinModel{Float64}:
  Frequency (Hz)      : 10.0
  DC                  : 1.0
  Sine amplitude (Q)  : -0.5
  Cosine amplitude (I): 1.2

julia> SinModel(1, 0, torect(1, -π/2)...)  # A pure sine
Sinusoidal parameters SinModel{Float64}:
  Frequency (Hz)      : 1.0
  DC                  : 0.0
  Sine amplitude (Q)  : 1.0
  Cosine amplitude (I): 0.0
```
"""
SinModel(f, DC, Q, I) = SinModel(promote(f, DC, Q, I)...)

"""
    SinModel{T}(; f, DC, Q, I)

Construct a `SinModel{T}` specifying each parameter by name.


Example
=======

```
julia> SinModel(DC = 1, Q = -0.5, I = 1.2, f = 10)
Sinusoidal parameters SinModel{Float64}:
  Frequency (Hz)      : 10.0
  DC                  : 1.0
  Sine amplitude (Q)  : -0.5
  Cosine amplitude (I): 1.2
```
"""
SinModel(; f, DC, Q, I) = SinModel(f, DC, Q, I)

"""
    (P::SinModel)(t)

Evaluate the sinusoidal function specified by the parameters `P` at the values given
by `t`, which may be a scalar or a collection.

Example
=======

```
julia> P = SinModel(DC = 1, Q = -0.5, I = 1.2, f = 10)
julia> t = range(0, 0.7, length = 5)
julia> P(t)
5-element Vector{Float64}:
  2.2
  1.4999999999999996
 -0.2000000000000004
  0.49999999999999944
  2.200000000000001
```
"""
(params::SinModel)(t) = @. params.Q*sin(2π*params.f*t) + params.I*cos(2π*params.f*t) + params.DC
#TODO: Remove implicit broadcasting

function Base.show(io::IO, params::SinModel{T}) where {T}
    println(io, "Sinusoidal parameters SinModel{$T}:")
    println(io, "  Frequency (Hz)      : $(params.f)")
    println(io, "  DC                  : $(params.DC)")
    println(io, "  Sine amplitude (Q)  : $(params.Q)")
    println(io, "  Cosine amplitude (I): $(params.I)")
end

"""
    MixedLinSinModel{T <: Real} <: SRModel

Parameters `f`, `DC`, `Q`, `I` and `m` of the sinusoidal function
``s(x) = DC + mx + I\\cos(2πfx) + Q\\sin(2πfx)``.

See also: [`SinModel`](@ref).
"""
struct MixedLinSinModel{T <: Real} <: SRModel
    f  :: T  # frequency in Hz
    DC :: T
    Q  :: T  # sine amplitude
    I  :: T  # cosine amplitude
    m  :: T  # linear term
end

"""
    MixedLinSinModel{T}(f, DC, Q, I, m) <: SRModel

Construct a `MixedLinSinModel{T}` with the given parameters, promoting to a common type
`T` if necessary.

Example
=======

```
julia> MixedLinSinModel(10, 0, -1.2, 0.4, 2.1)
Mixed Linear-Sinusoidal parameters MixedLinSinModel{Float64}:
  Frequency (Hz)       : 10.0
  DC                   : 0.0
  Sine amplitude (Q)   : -1.2
  Cosine amplitude (I) : 0.4
  Linear term (m)      : 2.1
```
"""
MixedLinSinModel(f, DC, Q, I, m) = MixedLinSinModel(promote(f, DC, Q, I, m)...)

"""
    MixedLinSinModel(; f, DC, Q, I, m)

Construct a `MixedLinSinModel{T}` specifying each parameter by name.

Example
=======

```
julia> MixedLinSinModel(f = 10, DC = 0, m = 2.1, Q= -1.2, I = 0.4)
Mixed Linear-Sinusoidal parameters MixedLinSinModel{Float64}:
  Frequency (Hz)       : 10.0
  DC                   : 0.0
  Sine amplitude (Q)   : -1.2
  Cosine amplitude (I) : 0.4
  Linear term (m)      : 2.1
```
"""
MixedLinSinModel(; f, DC, Q, I, m) = MixedLinSinModel(f, DC, Q, I, m)

"""
    (P::MixedLinSinModel)(t)

Evaluate the mixed linear-sinusoidal function specified by the parameters `P` at
the values given by `t`, which may be a scalar or a collection.

Example
=======

```
julia> P = MixedLinSinModel(f = 10, DC = 0, m = 2.1, Q= -1.2, I = 0.4)
julia> t = range(0, 0.7, length = 5)
julia> P(t)
5-element Vector{Float64}:
  0.4
  1.5674999999999997
  0.3349999999999989
 -0.09750000000000014
  1.870000000000002
```
"""
(params::MixedLinSinModel)(t) = @. params.Q*sin(2π*params.f*t) + params.I*cos(2π*params.f*t) +
                                       params.m*t + params.DC

function Base.show(io::IO, params::MixedLinSinModel{T}) where {T}
    println(io, "Mixed Linear-Sinusoidal parameters MixedLinSinModel{$T}:")
    println(io, "  Frequency (Hz)       : $(params.f)")
    println(io, "  DC                   : $(params.DC)")
    println(io, "  Sine amplitude (Q)   : $(params.Q)")
    println(io, "  Cosine amplitude (I) : $(params.I)")
    println(io, "  Linear term (m)      : $(params.m)")
end

### not yet implemented

"""
    GenSinProblem <: SinusoidalFunctionParameters

Not yet implemented.

See also: [`SinProblem`](@ref), [`MixedLinSinModel`](@ref), [`DampedSinProblem`](@ref)
(not yet implemented).
"""
struct GenSinProblem
end

"""
    DampedSinProblem <: SinusoidalFunctionParameters

Not yet implemented.

See also: [`SinProblem`](@ref), [`MixedLinSinModel`](@ref), [`GenSinProblem`](@ref)
(not yet implemented).
"""
struct DampedSinProblem
end
