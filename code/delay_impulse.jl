# Demonstration of high-accuracy simulation of delay systems
using ControlSystems
using OrdinaryDiffEq, DelayDiffEq, Statistics, LinearAlgebra, Polynomials
using Plots
f = (14, "Computer Modern")
gr(titlefont=f, tickfont=f,legendfont=f,guidefont=f,grid=false)

##
# Computes the time-response of feedback(1/s, delay(1.0)) up to numerical precision
function dde_example_exact(t)
    M = ceil(Int, t[end]) # Number of time segments
    
    # Compute matrix H where element (k, j) correponds
    # to the coefficient of (t-k)^(j-1) in the kth time segment
    #H = zeros(BigFloat, M, M)
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
        y[inds] .= pk.(t[inds] .- big(k-1))
    end
    
    return y'
end

t = 0:0.1:10
ytrue = dde_example_exact(big.(t))
plot(t, ytrue', label="\$ 1/(s+e^{-s}) \$")
savefig("impulse-example.pdf")

## Code part 1
s = tf("s")
sys = feedback(1/s, delay(1))
y, = impulse(sys, t)
maximum(abs,ytrue-y) # 8.9e-8

## Code part 2
abstol = reltol = 1e-12
y,     = impulse(sys, t; abstol, reltol)
maximum(abs, ytrue-y) # 4.5e-13
## END Code part 2
median((ytrue-y)./y)

## Code part 3
sys    = feedback(BigFloat(1)/s,delay(1))
abstol = reltol = 1e-30
alg    = MethodOfSteps(Vern9())
y,     = impulse(sys, big.(t); alg, abstol, reltol)
maximum(abs, ytrue-y) # 2.4e-32
## END Code part 3
median(max.(min.(abs.((ytrue-y)./y),BigFloat(1e40)),BigFloat(1e-40)))

