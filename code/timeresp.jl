using ControlSystems
s   = tf("s")
sys = feedback(1/s,delay(1))
t   = 0:0.1:20
impulseplot(sys, t)

using OrdinaryDiffEq, DelayDiffEq
alg    = MethodOfSteps(Tsit5())
abstol = reltol = 1e-12
y,     = impulse(sys, t;  alg, abstol, reltol)
maximum(abs, ytrue-y) # 2.23e-12
sys    = feedback(BigFloat(1)/s,delay(1))
abstol = reltol = 1e-25
alg    = MethodOfSteps(Vern6())
y,     = impulse(sys, t; alg, abstol, reltol)
maximum(abs, ytrue-y) # 1.53e-25