using ControlSystems
using GenericLinearAlgebra, OrdinaryDiffEq, DelayDiffEq, Statistics, LinearAlgebra, Polynomials
using Plots
f = (14, "Computer Modern")
gr(titlefont=f, tickfont=f,legendfont=f,guidefont=f,grid=false)


# Computes the time-response of feedback(1/s, delay(1.0)) up to numerical precision
function dde_example_exact(t)
  M = ceil(Int, t[end]) # Number of time segments

  # Compute matrix H where element (k, j) correponds
  # to the coefficient of (t-k)^(j-1) in the kth time segment
  H = zeros(BigFloat, M, M)

  H[1,1] = 1

  for k=2:M
      H[k, 1] = sum(H[k-1, :]) # Compute the constant (previous segment evaluated at 1)
      H[k, 2:end] = -H[k-1, 1:end-1] ./ (1:M-1) # Integrate the time-delayed feedback
  end

  # Evaluate the polynomial solution in each segment
  y = zeros(BigFloat, size(t))
  for k=1:M
    inds =  k-1 .<= t .<= k
    pk = Polynomial(H[k,:])
    y[inds] .= pk.(t[inds] .- (k-1))
  end

  return y'
end

t = 0:0.1:20
ytrue = dde_example_exact(t)
### Code part 1
s = tf("s")
sys = feedback(1/s,delay(1))
impulseplot(sys, t)
### END Code part 1
savefig("impulse-example.pdf")

y, = impulse(sys, t)

maximum(abs,ytrue-y) # 1e-7
median((ytrue-y)./y) # 1e-8

### Code part 2
alg    = MethodOfSteps(Tsit5())
abstol = reltol = 1e-12
y,     = impulse(sys, t;  alg, abstol, reltol)
maximum(abs, ytrue-y) # 1.75e-11
### END Code part 2
median((ytrue-y)./y)

### Code part 3
sys    = feedback(BigFloat(1)/s,delay(1))
abstol = reltol = 1e-25
alg    = MethodOfSteps(Vern6())
y,     = impulse(sys, t; alg, abstol, reltol)
maximum(abs, ytrue-y) # 2.6e-17
### END Code part 3
median(max.(min.(abs.((ytrue-y)./y),BigFloat(1e40)),BigFloat(1e-40)))
