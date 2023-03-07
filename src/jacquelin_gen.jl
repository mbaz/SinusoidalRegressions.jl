"""
    dampedsinfit_j :: DampedSinusoidP

Note: This function is not yet implemented.

Perform a damped sinusoid fit of the independent variables `X` and dependent
variables `Y`, using the method of integral equations.

The data is fit to the model

``f(x, f, Q, I, a) = exp(ax)(Qsin(2πf) + Icos(2πf))``

where ``a < 0``.  No initial guess as to the values of the parameters is needed. The values
in `X` must be sorted in ascending order.

The regression method implemented here is based on J. Jacquelin, "Régressions et équations
intégrales", 2014 (available at
https://fr.scribd.com/doc/14674814/Regressions-et-equations-integrales)
"""
function dampedsinfit_j(X, Y)
    println("Not yet implemented")
end
