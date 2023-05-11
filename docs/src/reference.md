# Reference

This package's workhorse function is `sinfit(p::SRProblem, a::SRAlgorithm)`. It takes a problem `p`
and calculates a fit using algorithm `a`.

```@docs
sinfit
```

## Algorithms

### IEEE 1057

Sinusoidal regressions specified by 1057-2017, IEEE Standard for Digitizing
Waveform Recorders. These are the same algorithms specified in 1241-2010, IEEE
Standard for Terminology and Test Methods for Analog-to-Digital Converters.

```@docs
IEEE1057
```

### Method of integral equations

Four types of sinusoidal regressions using the algorithms proposed by J.
Jacquelin in "Régressions et Equations Intégrales". These algorithms are not
iterative, and they do not require an initial estimate of the frequency or any
other parameter. They may be used by themselves, or to calculate initial
estimates for more-precise least-squares methods based on non-linear
optimization (described below).

```@docs
IntegralEquations
```

### Non-linear optimization

Wrapper for `LsqFit.curve_fit()`, which defines the appropriate model,
calculates initial parameter estimates if not provided by the user, and wraps the
returned fit in the appropriate parameter container.

```@docs
LevMar
```

### Liang

Algorithm designed for fitting when only a fraction of a period is sampled.

```@docs
Liang
```

## Problems

An `SRProblem` encapsulates what is known: at a minimum, the sampling times and the samples, and
possibly also estimates or bounds on some or all parameters. The problem type specifies the
number of unknown parameters, and the model (sinusoidal or mixed linear-sinusoidal).

```@docs
Sin3Problem
Sin4Problem
MixedLinSin4Problem
MixedLinSin5Problem
```

## Problem-Algorithm matrix

The following table shows which algorithms can solve which sinusoidal regression problems:

|                         |`IEEE1057` | `IntegralEquations` | `LevMar` | `Liang` |
|-------------------------|:---------:|:-------------------:|:--------:|:-------:|
|`Sin3Problem`            | ✓         |                     | ✓        |         |
|`Sin4Problem`            | ✓         | ✓                   | ✓        | ✓       |
|`MixedLinearSin4Problem` |           | ✓                   | ✓        |         |
|`MixedLinearSin5Problem` |           | ✓                   | ✓        |         |

## Sinusoidal parameters

```@docs
SinModel
SinModel(::Any, ::Any, ::Any, ::Any)
SinModel(; ::Any, ::Any, ::Any, ::Any)
SinModel(::Any)
MixedLinSinModel
MixedLinSinModel(::Any, ::Any, ::Any, ::Any, ::Any)
MixedLinSinModel(; ::Any, ::Any, ::Any, ::Any, ::Any)
MixedLinSinModel(::Any)
#SinusoidalRegressions.GenSinModel
#SinusoidalRegressions.DampedSinModel
```

## Error measurement

Functions to easily compare a calculated fit to actual data, using root mean-square errors,
or mean absolute errors.

```@docs
rmse
mae
```

## Plotting

This package provides recipes for `Plots.jl`. This makes it very easy to plot data and the
calculated fit.

There are two recipes included: plotting one fit, and plotting data and up to two fits.

### Plotting a fit

```julia
plot(X, fit::T, ; fitlabel = "") where {T <: SinusoidalFunctionParameters}
```

Plot the model with parameters given by `fit` at the points specified in collection `X`. The
optional keyword argument `fitlabel` controls the label (legend) of the plot.

##### Example

```@example 2
using SinusoidalRegressions
using Random: seed! # hide
using Plots # hide
seed!(5678) # hide
X = range(0, 1, length = 100)
Y = 2 .+ 3cos.(2*pi*5*X) .- 0.2sin.(2*pi*5*X) .+ 0.1*randn(100)
fit = sinfit(Sin3Problem(X, Y, 5), IEEE1057())
plot(X, fit, fitlabel = "example")
```

### Plotting data and one or two different fits

```julia
plot(X, Y, fit::T ; exact :: Union{T, Nothing} = nothing,
                    samples    = 100,
                    fitlabel   = "fit",
                    datalabel  = "data",
                    exactlabel = "exact") where {T <: SinusoidalFunctionParameters}
```
Plot the data `Y` and the model with parameters given by `fit` evaluated at the points specified
in the collection `X`. The following keyword arguments are supported (with default values in
parenthesis):
* `exact::T = nothing`: add one more model to the plot. This allows plotting the exact function, the data and the fit in a single plot.
* `samples::Int = 100`: `fit` and `exact` are evaluated for `range(first(X), last(X), length = samples)`.
* `fitlabel::String = "fit"`: the label (legend) for the fit plot.
* `datalabel::String = "data"`: the label for the data plot.
* `exactlabel::String = "exact"`: the label for the exact plot.

##### Example

```@example 2
X = range(0, 1, length = 20)
exact = SinModel(5, 2, -0.2, 3)
Y = exact.(X) .+ 0.3*randn(20)
fit = sinfit(Sin4Problem(X, Y), IntegralEquations())
plot(X, Y, fit, exact = exact,
     datalabel = "measurements",
     fitlabel = "NLS",
     exactlabel = "true")
```

## Abstract types

```@docs
SRModel
SRAlgorithm
SRProblem
```
