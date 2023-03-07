"""
    gensinfit_j :: GenSinusoidP

Note: This function is not yet implemented.

Perform a general sinusoidal fit of the independent variables `X` and dependent
variables `Y`, using the method of integral equations.

The data is fit to the model

``f(x, f, Q, I, L...) = Qsin(2πf) + Icos(2πf) + L[1]F[1](x) + ... + L[m]F[m](x)``

where `F` is a collection of known functions and `L` is a collection of scalars. No initial
guess as to the values of the parameters is needed. The values in `X` must be sorted in
ascending order.

The regression method implemented here is based on J. Jacquelin, "Régressions et équations
intégrales", 2014 (available at
https://fr.scribd.com/doc/14674814/Regressions-et-equations-integrales)
"""
function gensinfit_j(X, Y, F)
    println("Not yet implemented")
end
