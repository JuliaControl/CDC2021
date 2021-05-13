trig_identity(x) = sin(x)^2 + cos(x)^2
println(trig_identity(0.5))

using Symbolics
@variables x
trig_identity(x)
println(simplify(trig_identity(x)))

X = randn(2,2)   # Random 2×2 matrix
println(trig_identity(X))
