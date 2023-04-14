
function liang(prob::Sin4Problem, a::Liang)
    (; X, Y) = prob
    (;threshold, iterations, q) = a
    N = length(X)
    length(Y) ≠ N && error("Input arguments should have the same length (got $N and $(length(Y)))")

    # step 2
    τ = last(X) - first(X)  # observation interval length
    v = (N-1)/τ             # average sampling rate
    f0 = 1/τ                # true frequency is smaller than f0, since we assume P < 1
    while true
        if f0 > q/τ
            break
        else
            q = q/2
        end
    end

    # step 3
    fL = q/τ
    fR = 2/τ

    # step 4
    fM = fL + 0.618*(fR - fL)
    fT = fR - 0.618*(fR - fL)

    count = 0
    while true
        # step 5
        fitL = _ieee1057_3P(X, Y, fL)
        ρL = rmse(fitL, Y, X)
        fitR = _ieee1057_3P(X, Y, fR)
        ρR = rmse(fitR, Y, X)
        fitM = _ieee1057_3P(X, Y, fM)
        ρM = rmse(fitM, Y, X)
        fitT = _ieee1057_3P(X, Y, fT)
        ρT = rmse(fitT, Y, X)

        # step 6
        if ρM < ρT
            ρ = ρM
            fL = fT
            fT = fM
            fM = fL + 0.618*(fR - fL)
        else
            ρ = ρT
            fR = fM
            fM = fT
            fT = fR - 0.618*(fR - fL)
        end

        # step 7
        if abs((ρM - ρT)/ρT) < threshold
            if ρM < ρT
                return fitM
            else
                return fitT
            end
        end
        count += 1
        if count > iterations
            @warn "Maximum iterations reached" iterations
            if ρM < ρT
                return fitM
            else
                return fitT
            end
        end
    end
end
