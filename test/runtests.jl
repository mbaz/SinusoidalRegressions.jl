using SinusoidalRegressions
using Test

@testset "Example 1 from Jacquelin's paper: sinusoidal regression" begin
    # Example is on page 23; expected results are on page 32
    X = [-1.983, -1.948, -1.837, -1.827, -1.663, -0.815,
         -0.778, -0.754, -0.518,  0.322,  0.418,  0.781,
          0.931,  1.510,  1.607]
    Y = [ 0.936,  0.810,  0.716,  0.906,  0.247, -1.513,
         -1.901, -1.565, -1.896,  0.051,  0.021,  1.069,
          0.862,  0.183,  0.311]

    # Jacquelin's method
    (;f, DC, Q, I) = sinfit_j(X, Y)
    @test 2*pi*f ≈  2.02074  atol = 0.00001
    @test DC     ≈ -0.405617 atol = 0.00001
    @test Q      ≈  1.2752   atol = 0.0001
    @test I      ≈ -0.577491 atol = 0.000001

    # IEEE 1057 3-parameter
    (;DC, Q, I) = ieee1057(X, Y, 2.02074/(2*pi))
    @test DC ≈ -0.405617 atol = 0.00001
    @test Q  ≈  1.2752   atol = 0.0001
    @test I  ≈ -0.577491 atol = 0.00001

    # IEEE 4-parameter (algorithm fails to converge)
    @test SinusoidalRegressions.ieee1057_testconvergence(X, Y, f = 2.02074/(2*pi)) == false

    # LsqFit
    (;f, DC, Q, I) = sinfit(X, Y)
    @test 2*pi*f ≈ 2.02074   atol = 0.04
    @test DC     ≈ -0.405617 atol = 0.02
    @test Q      ≈  1.2752   atol = 0.02
    @test I      ≈ -0.577491 atol = 0.01
end

@testset "Polar notation for SinusoidP" begin
    s = SinusoidP(1, 0, 1, 0, polar = true) # A pure cosine
    @test s.f == 1
    @test s.DC == 0
    @test s.Q ≈ 0
    @test s.I ≈ 1
    s = SinusoidP(1, 0, 1, -π/2, polar = true) # A pure sine
    @test s.f == 1
    @test s.DC == 0
    @test s.Q ≈ 1 atol = 1e-12
    @test s.I ≈ 0 atol = 1e-12
    s = SinusoidP(1, 5, 1, -π/4, polar = true) # A mix
    @test s.f == 1
    @test s.DC == 5
    @test s.Q ≈ sqrt(2)/2 atol = 1e-12
    @test s.I ≈ sqrt(2)/2 atol = 1e-12
end

@testset "Easy sinusoidal regression tests" begin
    # Easy tests: lots of equidistant points, no noise.
    t = range(0, 1, length=100)

    # Q = 0
    exact = SinusoidP(f = 3.0, DC = 1.25, I = 4, Q = 0)
    Y = exact.(t)

    # integral equations
    (;f, DC, Q, I) = sinfit_j(t, Y)
    @test f  ≈ exact.f  atol = 0.0001
    @test DC ≈ exact.DC atol = 0.0001
    @test Q  ≈ exact.Q  atol = 0.0001
    @test I  ≈ exact.I  atol = 0.0001

    # ieee 1057 3-parameter
    (;DC, Q, I) = ieee1057(t, Y, exact.f)
    @test DC ≈ exact.DC atol = 0.0001
    @test Q  ≈ exact.Q  atol = 0.0001
    @test I  ≈ exact.I  atol = 0.0001

    # ieee 1057 4-parameter
    (;f, DC, Q, I) = ieee1057(t, Y)
    @test SinusoidalRegressions.ieee1057_testconvergence(t, Y) == true
    @test f  ≈ exact.f  atol = 0.001
    @test DC ≈ exact.DC atol = 0.001
    @test Q  ≈ exact.Q  atol = 0.001
    @test I  ≈ exact.I  atol = 0.001

    # LsqFit
    (;f, DC, Q, I) = sinfit(t, Y)
    @test f  ≈ exact.f  atol = 0.0001
    @test DC ≈ exact.DC atol = 0.0001
    @test Q  ≈ exact.Q  atol = 0.0001
    @test I  ≈ exact.I  atol = 0.0001

    # I = 0
    exact = SinusoidP(f = 3.0, DC = 1.25, I = 0, Q = -2)
    Y = exact.(t)

    # integral equations
    (;f, DC, Q, I) = sinfit_j(t, Y)
    @test f  ≈ exact.f  atol = 0.0001
    @test DC ≈ exact.DC atol = 0.0001
    @test Q  ≈ exact.Q  atol = 0.0001
    @test I  ≈ exact.I  atol = 0.0001

    # ieee 1057 3-parameter
    (;DC, Q, I) = ieee1057(t, Y, exact.f)
    @test DC ≈ exact.DC atol = 0.0001
    @test Q  ≈ exact.Q  atol = 0.0001
    @test I  ≈ exact.I  atol = 0.0001

    # ieee 1057 4-parameter
    (;f, DC, Q, I) = ieee1057(t, Y)
    @test SinusoidalRegressions.ieee1057_testconvergence(t, Y) == true
    @test f  ≈ exact.f  atol = 0.0001
    @test DC ≈ exact.DC atol = 0.0001
    @test Q  ≈ exact.Q  atol = 0.0001
    @test I  ≈ exact.I  atol = 0.0001

    # LsqFit
    (;f, DC, Q, I) = sinfit(t, Y)
    @test f  ≈ exact.f  atol = 0.0001
    @test DC ≈ exact.DC atol = 0.0001
    @test Q  ≈ exact.Q  atol = 0.0001
    @test I  ≈ exact.I  atol = 0.0001

    # I != 0, Q != 0  -- check less digits
    exact = SinusoidP(f = 2.0, DC = 0.25, I = 0.48, Q = -1.5)
    Y = exact.(t)

    # integral equations
    (;f, DC, Q, I) = sinfit_j(t, Y)
    @test f  ≈ exact.f  atol = 0.01
    @test DC ≈ exact.DC atol = 0.01
    @test Q  ≈ exact.Q  atol = 0.01
    @test I  ≈ exact.I  atol = 0.01

    # ieee 1057 3-parameter
    (;DC, Q, I) = ieee1057(t, Y, exact.f)
    @test DC ≈ exact.DC atol = 0.01
    @test Q  ≈ exact.Q  atol = 0.01
    @test I  ≈ exact.I  atol = 0.01

    # ieee 1057 4-parameter (fails in this case)
    @test SinusoidalRegressions.ieee1057_testconvergence(t, Y) == false

    # LsqFit -- more exact, so check more digits
    (;f, DC, Q, I) = sinfit(t, Y)
    @test f  ≈ exact.f  atol = 0.0001
    @test DC ≈ exact.DC atol = 0.0001
    @test Q  ≈ exact.Q  atol = 0.0001
    @test I  ≈ exact.I  atol = 0.0001
end

@testset "IEEE 1057 4-Parameter specific tests" begin
    t = range(0, 4, length = 1337)
    exact = SinusoidP(f = 2.0, DC = 0.25, I = 0.48, Q = -1.5)
    Y = exact.(t)
    @test SinusoidalRegressions.ieee1057_testconvergence(t, Y) ==  true
    @test SinusoidalRegressions.ieee1057_testconvergence(t, Y, f = 1.9) ==  false

    (;f, DC, Q, I) = ieee1057(t, Y)
    @test f  ≈ exact.f  atol = 0.0001
    @test DC ≈ exact.DC atol = 0.0001
    @test Q  ≈ exact.Q  atol = 0.0001
    @test I  ≈ exact.I  atol = 0.0001

    t = range(0, 1, length = 2313)
    exact = SinusoidP(f = 3.0, DC = -0.98, I = 2.15, Q = 3.87)
    Y = exact.(t)
    (;f, DC, Q, I) = ieee1057(t, Y)
    @test f  ≈ exact.f  atol = 0.0001
    @test DC ≈ exact.DC atol = 0.0001
    @test Q  ≈ exact.Q  atol = 0.0001
    @test I  ≈ exact.I  atol = 0.0001
end

@testset "Example 2 from Jacquelin's paper: mixed linear-sinusoidal regression" begin
    # First, "large scatter" example from Jacquelin's paper -- page 51
    X = [0, 0.31, 0.6, 0.71, 0.82, 1.62, 2.03, 2.73, 2.93, 3.19,
         3.28, 3.68, 3.72, 4.26, 4.75, 6.72, 7.75, 8.41, 8.61, 9.17]
    Y = [1.35, 1.95, 2.20, 1.91, 1.88, 2.56, 2.74, 2.33, 2.82, 2.14,
         2.29, 1.76, 2.46, 1.70, 2.20, 4.43, 5.59, 7.03, 7.19, 6.84]
    actual = MixedLinearSinusoidP(0.8/(2π), 0.4, 1.2, 0.75, 0.6)
    expected = MixedLinearSinusoidP(0.861918/(2π), 0.571315, 1.113310, 0.613981, 0.582827)
    (;f, DC, Q, I, m) = mixlinsinfit_j(X, Y)
    @test f  ≈ expected.f  atol = 0.01
    @test DC ≈ expected.DC atol = 0.01
    @test Q  ≈ expected.Q  atol = 0.02
    @test I  ≈ expected.I  atol = 0.02
    @test m  ≈ expected.m  atol = 0.01

    # LsqFit
    (;f, DC, Q, I, m) = mixlinsinfit(X, Y)
    @test f  ≈ actual.f  atol = 0.005
    @test DC ≈ actual.DC atol = 0.15
    @test Q  ≈ actual.Q  atol = 0.2
    @test I  ≈ actual.I  atol = 0.01
    @test m  ≈ actual.m  atol = 0.01

    # Second, "small scatter" example from Jacquelin's paper -- page 53
    X = range(start = 0.5, step = 0.5, length = 20)
    Y = [1.66, 2.20, 2.83, 2.66, 2.59, 2.53, 2.12, 1.85, 1.85, 1.97,
         2.16, 2.86, 3.42, 4.56, 5.11, 6.00, 6.91, 7.16, 7.56, 7.41]
    actual = MixedLinearSinusoidP(0.8/(2π), 0.4, 1.2, 0.75, 0.6)
    expected = MixedLinearSinusoidP(0.809132/(2π), 0.266603, 1.278916, 0.680126, 0.613464)
    (;f, DC, Q, I, m) = mixlinsinfit_j(X, Y)
    @test f  ≈ expected.f  atol = 0.001
    @test DC ≈ expected.DC atol = 0.001
    @test Q  ≈ expected.Q  atol = 0.001
    @test I  ≈ expected.I  atol = 0.001
    @test m  ≈ expected.m  atol = 0.001

    # LsqFit
    (;f, DC, Q, I, m) = mixlinsinfit(X, Y)
    @test f  ≈ actual.f  atol = 0.002
    @test DC ≈ actual.DC atol = 0.15
    @test Q  ≈ actual.Q  atol = 0.1
    @test I  ≈ actual.I  atol = 0.1
    @test m  ≈ actual.m  atol = 0.015
end

@testset "RMSE tests" begin

end
