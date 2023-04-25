# Introduction

SinusoidalRegressions.jl aims to provide a set of functions for conveniently
fitting noisy data to a variety of sinusoidal models, with or without initial
estimates of the parameters.

## Features

This package supports five sinusoidal models:

* 3-parameter sinusoidal, where the frequency ``f`` is known exactly:

  ``s_s(x; \textrm{DC}, I, Q) = \textrm{DC} + Q\sin(2πfx) + I\cos(2πfx)``

* 4-parameter sinusoidal, where the frequency ``f`` is unknown:

  ``s_s(x; f, \textrm{DC}, I, Q) = \textrm{DC} + Q\sin(2πfx) + I\cos(2πfx)``

* Mixed linear-sinusoidal:

  ``s_{mls}(x; f, \textrm{DC}, I, Q, m) = \textrm{DC} + mx + Q\sin(2πfx) + I\cos(2πfx)``

* (Not yet implemented) General sinusoidal with unknown coefficients ``λ_1, λ_2, \ldots, λ_m``, and known functions ``g_1(x), g_2(x), \ldots, g_m(x)``:

  ``s_g(x; f, I, Q, λ_n) = Q\sin(2πfx) + I\cos(2πfx) + λ_1g_1(x) + \ldots + λ_mg_m(x)``

* (Not yet implemented) Damped sinusoidal:

  ``s_d(x; f, I, Q, α) = \exp(αx) ( Q\sin(2πfx) + I\cos(2πfx) )``

!!! note "A note on the frequency"
    The frequency ``f`` is assumed to be expressed in hertz throughout.

Four types of sinusoidal fitting algorithms are provided:

  * IEEE 1057 3-parameter and 4-parameter algorithms[^1]. These are able to fit only the sinusoidal model
    ``s_s(x)``. The 4-parameter version requires a very close estimate of the frequency ``f`` and
    dense sampling of several periods of the sinusoid.

  * The algorithms proposed by J. Jacquelin, based on the idea of finding an integral
    equation whose solution is the desired sinusoidal model[^2]. This integral equation is linear
    and can be solved using linear least-squares. These algorithms have several benefits:

      * They are non-iterative.
      * They do not require initial estimates of any of the model's parameters.
      * They often perform quite well, but may also be used to calculate a good initial guess for the
        more powerful non-linear methods described below.

  * The non-linear fitting function `curve_fit` from the package
    [`LsqFit`](https://github.com/JuliaNLSolvers/LsqFit.jl). This function uses the Levenberg-Marquardt
    algorithm[^3][^4][^5] and is quite powerful, but it requires an initial estimate of the parameters.
    SinusoidalRegressions package "wraps" `curve_fit`, automatically defining the correct model
    and calculating initial parameter estimates (using IEEE 1057 or Jacquelin's algorithms, as
    appropriate) if none are provided by the user.

   * The algorithm proposed by Liang et al[^6], designed for the case where only a fraction of a period
     of a sinusoid is sampled.

[^1]: IEEE, "IEEE Standard for Digitizing Waveform Recorders", 2018.

[^2]:Jacquelin J., Régressions et Equations Intégrales, 2014, (online)
https://fr.scribd.com/doc/14674814/Regressions-et-equations-integrales.

[^3]: Levenberg, K., "A Method for the Solution of Certain Non-Linear Problems in Least Squares",
Quarterly of Applied Mathematics, 2 (2), 1944. doi:10.1090/qam/10666.

[^4]: Marquardt, D., "An Algorithm for Least-Squares Estimation of Nonlinear Parameters",
SIAM Journal on Applied Mathematics, 11 (2), 1944. doi:10.1137/0111030.

[^5] `LsqFit.jl` User Manual, (online) https://julianlsolvers.github.io/LsqFit.jl/latest/

[^6]: Liang et al, "Fitting Algorithm of Sine Wave with Partial Period Waveforms and
Non-Uniform Sampling Based on Least-Square Method." Journal of Physics:
Conference Series 1149.1 (2018)ProQuest. Web. 17 Apr. 2023.

## Installation

This package is in Julia's general registry. It can be installed by running
```julia
using Pkg; Pkg.add("SinusoidalRegressions)
```

## Roadmap and Contributions

The following features are on the roadmap towards version 1.0:

* Implement Jacquelin's algorithms for damped and general sinusoids.
* Enhance the front-end to `curve_fit` by allowing some paramters to be declared as known a priori, and
   fitting only the remaining parameters.
* Implement other algorithms and add support for other models, as suggested by the community.
* Improve tests.

Contributions, feature requests, bug reports and suggestions are welcome;
please use the [issue tracker on github](https://github.com/mbaz), or open a
discussion on the [Julia forums](https://discourse.julialang.org/).

## Tutorial

### Fitting experimental data

Fitting data using this package requires three steps:

1. Define a `Problem`. The problem encapsulates all known data: the sampling points and
   data at a minimum, and possibly also the frequency. Initial estimates and bounds are
   also part of the problem definition
2. Specify and optionally configure an `Algorithm`.
3. Run the `sinfit` function on the problem with the chosen algorithm. This function
   returns a `SinusoidalFunctionParameters` with the fitted parameters.

After the data is fit, it may be easily plotted, and its RMSE and MAE may be calculated.

As a first exmaple, let us fit noisy data to a mixed linear-sinusoidal model using Jacquelin's
integral equation algorithm. The model is

``s(x ; f, \textrm{DC}, Q, I, m) = \textrm{DC} + mx + Q\sin(2πfx) + I\cos(2πfx)``

with parameters ``f`` (the frequency is unknown), ``\textrm{DC}``, ``Q``, ``I`` and ``m``.

Our tasks are: define an "exact" model and generate noisy samples from it, fit the noisy samples,
determine the error, and plot the fit. First we need to learn about defining a model.

#### Model types

This package includes several types used to store the parameters of the different supported models:
* [`SinusoidP`](@ref) for representing sinusoidal models ``s_s(x)``.
* [`MixedLinearSinusoidP`](@ref) for representing mixed sinusoidal models ``s_{mls}(x)``.
* [`SinusoidalRegressions.GenSinusoidP`](@ref) for representing general sinusoidal models ``s_g(x)`` (not yet implemented).
* [`SinusoidalRegressions.DampedSinusoidP`](@ref) for representing damped sinusoidal models ``s_d(x)``.
All of these types are subtypes of the abstract type [`SinusoidalFunctionParameters`](@ref). 

Instances of these types are also function types (functors), making it easy to evaluate the
model at any point (or collection of points). A plot recipe is also provided. We will
use both of these features below.

Since we're assuming a mixed linear-sinusoidal model, we'll use the type `MixedLinearSinusoidP`.

#### Generating the "exact" data

We will first generate the "exact" or "true" data, and then add noise. Then, we'll fit the noisy
data and compare the fit to the exact model.

The exact parameters are ``f = 4.1``, ``\textrm{DC} = 0.2``, ``m = 1``, ``Q =
-0.75``, and ``I = 0.8``. We define the model as follows:
```@example 1
using SinusoidalRegressions
p_exact = MixedLinearSinusoidP(f = 4.1, DC = 0.2, m = 1, Q = -0.75, I = 0.8)
```
We now genereate 40 "exact" data points with a sampling rate of 40 Hz.
```@example 1
x = range(0, length = 40, step = 1/40)
d_exact = p_exact(x)  # exact data
nothing # hide
```
And finally, we add gaussian noise to the data:
```@example 1
using Random: seed!
seed!(6778899)

data = d_exact .+ 0.2*randn(40)  # noisy data
nothing # hide
```

#### Defining the problem



#### Fitting and error measurement

The vector `data` is assumed to contain the experimental or measured data, which
we wish to fit to the mixed linear-sinusoidal model. We do so as follows:
```@example 1
fit = mixlinsinfit_j(x, data)
```
We calculate the fit's root mean-square error and mean absolute error:
```@example 1
rmse(fit, p_exact, x)
```
```@example 1
mae(fit, p_exact, x)
```

#### Plotting

Plotting the measured data along with the fit:
```@example 1
using Plots
plot(x, data, fit)
```

#### Improving the fit with non-linear least squares

We may try to improve the fit by using `LsqFit.curve_fit` with the parameters estimated above. The
function `mixlinsinfit` calls `LsqFit.curve_fit` behind the scenes, setting up the model and calculating
the initial parameter estimates (using `mixlinsinfit_j`) if none are provided.
```@example 1
fit_nls = mixlinsinfit(x, data)
```
```@example 1
rmse(fit_nls, p_exact, x)
```
```@example 1
mae(fit_nls, p_exact, x)
```
As expected, this second fit has a smaller error. Plotting the data, the fit, and the exact
function for comparison:
```@example 1
plot(x, data, fit_nls, exact = p_exact)
```

### A more difficult scenario

Let us fit the same model as above, with fewer data points, more noise and non-equidistant samples.

```@example 1
x = sort(rand(20))
d_exact = p_exact(x)
data = d_exact .+ 0.3*randn(20)
fit = mixlinsinfit_j(x, data)    # using Jacquelin's algorithm
```
Error for Jacquelin's algorithm:
```@example 1
rmse(fit, p_exact, x)
```
```@example 1
mae(fit, p_exact, x)
```
As expected, we see larger errors than in the more benign scenario above. Let us improve the fit
using non-linear least-squares. Here, `fit` (calculated above) is explicitly given as an initial
guess:
```@example 1
fit_nls = mixlinsinfit(x, data, fit)
```
Error for `LsqFit.curve_fit`:
```@example 1
rmse(fit_nls, p_exact, x)
```
```@example 1
mae(fit_nls, p_exact, x)
```
We see that `LsqFit.curve_fit` was able to improve the fit significantly. Plotting the data,
the fit and the exact function:
```@example 1
plot(x, data, fit_nls, exact = p_exact)
```

## References

