## SinusoidalRegressions.jl

[![CI](https://github.com/mbaz/SinusoidalRegressions.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/mbaz/SinusoidalRegressions.jl/actions/workflows/ci.yml)

SinusoidalRegressions.jl aims to provide a set of functions for conveniently fitting noisy data to a variety of sinusoidal models, with or without initial estimates of the parameters. The package is quite usable in its current state, but is still in development. Support for more sinusoidal models will be added in the future, and API changes cannot be ruled out.

Its documentation is found [here](https://mbaz.github.io/SinusoidalRegressions.jl/stable/).

### Package features:

* An implementation of IEEE 1057 fitting algorithms for 3 and 4 parameters.
* An implementation of the fitting algorithms developed by J. Jacquelin, based on integral equations that can be solved numerically and whose solution provide the desired fit. These algorithms do not require an initial parameter estimate.
* A front-end to the non-linear fitting function `curve_fit` from the package [`LsqFit`](https://github.com/JuliaNLSolvers/LsqFit.jl). This function uses the Levenberg-Marquardt algorithm and is quite powerful, but it requires an initial estimate of the parameters.
* Support for sinusoidal and mixed linear-sinusoidal models.

In addition, the package provides functions to calculate the RMSE and MAE when the exact parameters are known, and plot recipes for convenient plotting.
