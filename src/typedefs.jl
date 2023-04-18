
#
# Algorithms
#

abstract type Algorithm end

"""
    IEEE1057([iterations = 6]) <: Algorithm

Define an instance of the IEEE 1057 sinusoidal fitting algorithm.

Optional argument `iterations` specifies how many iterations to run before the algorithm
stops. The default value is 6, which is the value recommended by the standard. This value
is only used when calculating a 4-parameter fit.
"""
Base.@kwdef struct IEEE1057 <: Algorithm
    iterations::Int = 6
end

"""
    IntegralEquations() <: Algorithm

Define an instance of the integral-equations sinusoidal fitting algorithm described by
J.  Jacquelin in "Régressions et équations intégrales", 2014 (available at
https://fr.scribd.com/doc/14674814/Regressions-et-equations-integrales).

This algorithm does not accept any configuration parameters.
"""
struct IntegralEquations <: Algorithm end

"""
    LevMar([use_ga = false]) <: Algorithm

Define an instance of the Levenberg-Marquardt sinusoidal fitting algorithm.

If the optional argument `use_ga` is set to `true`, the algorithm will use geodesic
acceleration to potentially improve its performance and accuracy. See
[`curve_fit`](https://github.com/JuliaNLSolvers/LsqFit.jl) for more details.
"""
Base.@kwdef struct LevMar <: Algorithm
    use_ga :: Bool = false
end

"""
    Liang([threshold = 0.15, iterations = 100, q = 1e-5]) <: Algorithm

Define an instance of the sinusoidal fitting algorithm described in Liang et al,
"Fitting Algorithm of Sine Wave with Partial Period Waveforms and Non-Uniform Sampling
Based on Least-Square Method." Journal of Physics: Conference Series 1149.1
(2018)ProQuest. Web. 17 Apr. 2023.

This algorithm is designed for scenarios where only a fraction of a period of the sinusoid
has been sampled. Its optional parameters `threshold` and `q` are described in the paper.
Additionally, a maximum number of iterations may be specified in `iterations`.
"""
Base.@kwdef struct Liang <: Algorithm
    threshold  :: Float64 = 0.15
    iterations :: Int     = 100
    q          :: Float64 = 1e-5
end

#
# Problems
#

abstract type Problem end

"""
    Sin3Problem(X, Y, f , [DC, Q, I, lb, ub]) <: Problem

Define a three-parameter sinusoidal regression problem.

The data is fit to the model ``f(x; DC, Q, I) = DC + Q\\sin(2πfx) + I\\cos(2πfx)``. The sampling
instants are given by `X`, and the samples by `Y`. The frequency `f` (in Hz) is assumed to
be known exactly.

When using a fitting algorithm that accepts initial parameter estimates, these may be given
by the optional keyword arguments `DC`, `Q` and `I` (default `missing`). Lower and upper bounds
may be specified in `lb` and `ub`, which must be vectors of length 3.

See also: [`Sin4Problem`](@ref)
"""
Base.@kwdef struct Sin3Problem{T1, T2, T3} <: Problem where
                                              {T1 <: AbstractVector, T2 <: AbstractVector, T3 <: Real}
    X  :: T1  # sampling times
    Y  :: T2  # sample values
    f  :: T3  # exact known frequency
    DC = missing  # initial estimates (may be missing)
    Q  = missing
    I  = missing
    lb = missing  # lower bounds
    ub = missing  # upper bounds
end

Sin3Problem(X, Y, f ; kwargs...) = Sin3Problem(; X, Y, f, kwargs...)

"""
    Sin4Problem(X, Y ; [f, DC, Q, I, lb, ub]) <: Problem

Define a four-parameter sinusoidal regression problem.

The data is fit to the model ``f(x; f, DC, Q, I) = DC + Q\\sin(2πfx) + I\\cos(2πfx)``.
The sampling instants are given by `X`, and the samples by `Y`.

When using a fitting algorithm that accepts initial parameter estimates, these may be given
by the optional keyword arguments `f`, `DC`, `Q` and `I` (default `missing`). Lower and upper
bounds may be specified in `lb` and `ub`, which must be vectors of length 4.

See also: [`Sin3Problem`](@ref)
"""
Base.@kwdef struct Sin4Problem{T1, T2} <: Problem where
                                          {T1 <: AbstractVector, T2 <: AbstractVector}
    X  :: T1  # sampling times
    Y  :: T2  # sample values
    f  = missing  # initial estimates (may be missing)
    DC = missing
    Q  = missing
    I  = missing
    lb = missing  # lower bounds
    ub = missing  # uppper bounds
end

Sin4Problem(X, Y ; kwargs...) = Sin4Problem(; X, Y, kwargs...)

"""
    MixedLinSin4Problem(X, Y, f ; [DC, Q, I, m, lb, ub])

Define a four-parameter mixed linear-sinusoidal regression problem.

The data is fit to the model ``f(x; DC, Q, I, m) = DC + Q\\sin(2πfx) + I\\cos(2πfx) + mx``.
The sampling instants are given by `X`, and the samples by `Y`. The frequency `f` (in Hz) is
assumed to be known exactly.

When using a fitting algorithm that accepts initial parameter estimates, these may be given
by the optional keyword arguments `DC`, `Q` `I` and `m` (default `missing`). Lower and upper
bounds may be specified in `lb` and `ub`, which must be vectors of length 4.

See also: [`MixedLinSin5Problem`](@ref)
"""
Base.@kwdef struct MixedLinSin4Problem{T1, T2, T3} <: Problem where
                                                      {T1 <: AbstractVector, T2 <: AbstractVector, T3 <: Real}
    X :: T1
    Y :: T2
    f :: T3
    DC = missing
    Q  = missing
    I  = missing
    m  = missing
    lb = missing
    ub = missing
end

MixedLinSin4Problem(X, Y, f ; kwargs...) = MixedLinSin4Problem(; X, Y, f, kwargs...)

"""
    MixedLinSin5Problem(X, Y ; [f, DC, Q, I, m, lb, ub])

Define a five-parameter mixed linear-sinusoidal regression problem.

The data is fit to the model ``f(x ; f, DC, Q, I, m) = DC + Q\\sin(2πfx) + I\\cos(2πfx) + mx``.
The sampling instants are given by `X`, and the samples by `Y`.

When using a fitting algorithm that accepts initial parameter estimates, these may be given
by the optional keyword arguments `f`, `DC`, `Q` `I` and `m` (default `missing`). Lower and upper
bounds may be specified in `lb` and `ub`, which must be vectors of length 5.

See also: [`MixedLinSin4Problem`](@ref)
"""
Base.@kwdef struct MixedLinSin5Problem{T1, T2} <: Problem where
                                                  {T1 <: AbstractVector, T2 <: AbstractVector, T3 <: Real}
    X :: T1
    Y :: T2
    f = missing
    DC = missing
    Q  = missing
    I  = missing
    m  = missing
    lb = missing
    ub = missing
end

MixedLinSin5Problem(X, Y ; kwargs...) = MixedLinSin5Problem(; X, Y, kwargs...)

"""
    SinusoidalFunctionParameters

Supertype for the parameters of the different kinds of sinusoidal models
supported by SinusoidalRegressions.

See also: [`SinusoidP`](@ref), [`MixedLinearSinusoidP`](@ref)
"""
abstract type SinusoidalFunctionParameters end

"""
    SinusoidP{T <: Real} <: SinusoidalFunctionParameters

Parameters `f`, `DC`, `Q` and `I` of the sinusoidal function
``s(x) = DC + I\\cos(2πfx) + Q\\sin(2πfx)``.

See also: [`MixedLinearSinusoidP`](@ref), [`GenSinusoidP`](@ref) (not yet implemented),
[`DampedSinusoidP`](@ref) (not yet implemented).
"""
struct SinusoidP{T <: Real} <: SinusoidalFunctionParameters
    f  :: T  # frequency in Hz
    DC :: T
    Q  :: T  # sine amplitude
    I  :: T  # cosine amplitude
end

"""
    SinusoidP{T}(f, DC, Q, I)

Construct a `SinusoidP{T}` with the given parameters, promoting to a common type `T` if
necessary.

To express the model as ``s(x) = M\\cos(2πfx + θ)``, use the `torect` function.

Examples
========

```
julia> SinusoidP(10, 1, -0.5, 1.2)
Sinusoidal parameters SinusoidP{Float64}:
  Frequency (Hz)      : 10.0
  DC                  : 1.0
  Sine amplitude (Q)  : -0.5
  Cosine amplitude (I): 1.2

julia> SinusoidP(1, 0, torect(1, -π/2)...)  # A pure sine
Sinusoidal parameters SinusoidP{Float64}:
  Frequency (Hz)      : 1.0
  DC                  : 0.0
  Sine amplitude (Q)  : 1.0
  Cosine amplitude (I): 0.0
```
"""
SinusoidP(f, DC, Q, I) = SinusoidP(promote(f, DC, Q, I)...)

"""
    SinusoidP{T}(; f, DC, Q, I)

Construct a `SinusoidP{T}` specifying each parameter by name, promoting to a common type `T` if
necessary.


Example
=======

```
julia> SinusoidP(DC = 1, Q = -0.5, I = 1.2, f = 10)
Sinusoidal parameters SinusoidP{Float64}:
  Frequency (Hz)      : 10.0
  DC                  : 1.0
  Sine amplitude (Q)  : -0.5
  Cosine amplitude (I): 1.2
```
"""
SinusoidP(; f, DC, Q, I) = SinusoidP(f, DC, Q, I)

"""
    (P::SinusoidP)(t)

Evaluate the sinusoidal function specified by the parameters `P` at the values given
by `t`, which may be a scalar or a collection.

Example
=======

```
julia> P = SinusoidP(DC = 1, Q = -0.5, I = 1.2, f = 10)
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
(params::SinusoidP)(t) = @. params.Q*sin(2π*params.f*t) + params.I*cos(2π*params.f*t) + params.DC

function Base.show(io::IO, params::SinusoidP{T}) where {T}
    println(io, "Sinusoidal parameters SinusoidP{$T}:")
    println(io, "  Frequency (Hz)      : $(params.f)")
    println(io, "  DC                  : $(params.DC)")
    println(io, "  Sine amplitude (Q)  : $(params.Q)")
    println(io, "  Cosine amplitude (I): $(params.I)")
end

"""
    MixedLinearSinusoidP{T <: Real} <: SinusoidalFunctionParameters

Parameters `f`, `DC`, `Q`, `I` and `m` of the sinusoidal function
``s(x) = DC + mx + I\\cos(2πfx) + Q\\sin(2πfx)``.

See also: [`SinusoidP`](@ref), [`GenSinusoidP`](@ref) (not yet implemented),
[`DampedSinusoidP`](@ref) (not yet implemented).
"""
struct MixedLinearSinusoidP{T <: Real} <: SinusoidalFunctionParameters
    f  :: T  # frequency in Hz
    DC :: T
    Q  :: T  # sine amplitude
    I  :: T  # cosine amplitude
    m  :: T  # linear term
end

"""
    MixedLinearSinusoidP{T}(f, DC, Q, I, m)

Construct a `MixedLinearSinusoidP{T}` with the given parameters, promoting to a common type
`T` if necessary.

Example
=======

```
julia> MixedLinearSinusoidP(10, 0, -1.2, 0.4, 2.1)
Mixed Linear-Sinusoidal parameters MixedLinearSinusoidP{Float64}:
  Frequency (Hz)       : 10.0
  DC                   : 0.0
  Sine amplitude (Q)   : -1.2
  Cosine amplitude (I) : 0.4
  Linear term (m)      : 2.1
```
"""
MixedLinearSinusoidP(f, DC, Q, I, m) = MixedLinearSinusoidP(promote(f, DC, Q, I, m)...)

"""
    MixedLinearSinusoidP(; f, DC, Q, I, m)

Construct a `MixedLinearSinusoidP{T}` specifying each parameter by name.

Example
=======

```
julia> MixedLinearSinusoidP(f = 10, DC = 0, m = 2.1, Q= -1.2, I = 0.4)
Mixed Linear-Sinusoidal parameters MixedLinearSinusoidP{Float64}:
  Frequency (Hz)       : 10.0
  DC                   : 0.0
  Sine amplitude (Q)   : -1.2
  Cosine amplitude (I) : 0.4
  Linear term (m)      : 2.1
```
"""
MixedLinearSinusoidP(; f, DC, Q, I, m) = MixedLinearSinusoidP(f, DC, Q, I, m)

"""
    (P::MixedLinearSinusoidP)(t)

Evaluate the mixed linear-sinusoidal function specified by the parameters `P` at
the values given by `t`, which may be a scalar or a collection.

Example
=======

```
julia> P = MixedLinearSinusoidP(f = 10, DC = 0, m = 2.1, Q= -1.2, I = 0.4)
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
(params::MixedLinearSinusoidP)(t) = @. params.Q*sin(2π*params.f*t) + params.I*cos(2π*params.f*t) +
                                       params.m*t + params.DC

function Base.show(io::IO, params::MixedLinearSinusoidP{T}) where {T}
    println(io, "Mixed Linear-Sinusoidal parameters MixedLinearSinusoidP{$T}:")
    println(io, "  Frequency (Hz)       : $(params.f)")
    println(io, "  DC                   : $(params.DC)")
    println(io, "  Sine amplitude (Q)   : $(params.Q)")
    println(io, "  Cosine amplitude (I) : $(params.I)")
    println(io, "  Linear term (m)      : $(params.m)")
end

### not yet implemented

"""
    GenSinusoidP <: SinusoidalFunctionParameters

Not yet implemented.

See also: [`SinusoidP`](@ref), [`MixedLinearSinusoidP`](@ref), [`DampedSinusoidP`](@ref)
(not yet implemented).
"""
struct GenSinusoidP
end

"""
    DampedSinusoidP <: SinusoidalFunctionParameters

Not yet implemented.

See also: [`SinusoidP`](@ref), [`MixedLinearSinusoidP`](@ref), [`GenSinusoidP`](@ref)
(not yet implemented).
"""
struct DampedSinusoidP
end
