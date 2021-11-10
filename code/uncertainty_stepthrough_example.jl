#=
In this example we will work through some typical usage of the toolbox. We will use an example of controller optimization with a plant model containing uncertain parameters.

We start by defining a nominal model, visalize it's Bode diagram and perform a time-domain simulation
=#
using ControlSystems, Plots
p = 1  
ζ = 0.3
ω = 1  
P = ss(tf([p*ω], [1, 2ζ*ω, ω^2]))
Ω = exp10.(-2:0.04:3)
kp,ki,kd =  1, 0.1, 0.1 # controller parameters

C  = ss(pid(; kp, ki, kd)*tf(1, [0.05, 1])^2) # Construct a PID controller with filter
G  = feedback(P*C) # Closed-loop system
S  = 1/(1 + P*C)   # Sensitivity function
CS = C*S           # Noise amplification
Gd = c2d(G, 0.1)   # Discretize the system
res = step(Gd,15)  # Step-response

mag = bode(S, Ω)[1][:]
plot(res, title="Time response", layout = 2, legend=:bottomright)
plot!(Ω, mag, title="Sensitivity function", xscale=:log10, yscale=:log10, subplot=2, legend=:bottomright, ylims=(3e-2, Inf))


#=
Next, we try to tune the controller using optimization. We therefore package the creation of the systems above in a function that takes in controller parameters and outputs a cost 
=#
using Optim, Statistics

const Msc = 1.2 # Constraint on Ms

function systems(params)
    kp,ki,kd = exp.(params)
    C  = ss(pid(; kp, ki, kd)*tf(1, [0.05, 1])^2)
    G  = feedback(P*C) # Closed-loop system
    S  = 1/(1 + P*C)   # Sensitivity function
    CS = C*S           # Noise amplification
    Gd = c2d(G, 0.1)   # Discretize the system
    res = step(Gd,15)  # Step-response
    C, G, S, CS, res
end

function cost(params)
    C, G, S, CS, res = systems(params)
    Ms = maximum(bode(S, Ω)[1]) # max sensitivity
    perf = mean(abs, 1 .- res.y)
    robust = (Ms > Msc ? 10000(Ms-Msc) : zero(eltype(params)))    
    noise = sum(bode(CS, Ω[end-30:end])[1])
    100perf + robust + 0.002noise
end

params  = log.([1,0.1,0.1]) # Initial guess
res = optimize(cost, params, Optim.Options(
    show_trace        = true,
    show_every        = 50,
))

## We can now perform the same computations as above to visualize the found controller
using Plots.Measures
𝕗 = (14, "Computer Modern")
gr(titlefont=𝕗, tickfont=𝕗, legendfont=𝕗, guidefont=𝕗, grid=false)
function plot_optimized(params, res)
    fig = plot(layout=2, size=(1000,400), bottommargin=2mm)
    for (i,params) = enumerate((params, res.minimizer))
        C, G, S, CS, r = systems(params)
        mag   = bode(S, Ω)[1][:]
        plot!(r, title="Time response", subplot=1, lab= i==1 ? "Initial" : "Optimized", legend=:bottomright)
        plot!(Ω, mag, title="Sensitivity function", xscale=:log10, yscale=:log10, subplot=2, lab= i==1 ? "Initial" : "Optimized", legend=:bottomright)
    end
    hline!([Msc], l=(:black, :dash), subplot=2, lab="", ylims=(3e-2, Inf))
    display(fig)
end

plot_optimized(params, res)

#=
Not bad :) 

The next step is to add uncertainty to the system. Lets say all the parameters (p, ζ, ω) are associated with a Gaussian uncertainty. We can create such parameters using MonteCarloMeasurements.jl
=#
using MonteCarloMeasurements: Particles, pmean, pmaximum, unsafe_comparisons
unsafe_comparisons(true)
±(m,s) = m + s*Particles(200)
p = 1   ± 0.1
ζ = 0.3 ± 0.05
ω = 1   ± 0.05
P = ss(tf([p*ω], [1, 2ζ*ω, ω^2]))

## If we now visualize the found controller, we find that it might violate the constraint we placed on Ms
plot_optimized(params, res)

#=
To alleviate this, we perform the optimization again, this time with the uncertain system.
In order to do this, we must change our cost function slightly so that it outputs a scalar instead of an uncertain numer
=#
cost(params)
#=
We do this by specifying which aspect of the probability distribution we are penalizing.

We add a call to pmaximum to say that we are interested in the worst-case Ms
and add a call to pmean to specify that we care about the mean performance
=#
function cost(params)
    C, G, S, CS, res = systems(params)
    Ms = maximum(bode(S, Ω)[1])     |> pmaximum# max sensitivity
    perf = mean(abs, 1 .- res.y)    |> pmean
    robust = (Ms > Msc ? 10000(Ms-Msc) : zero(eltype(params)))    
    noise = sum(bode(CS, Ω[end-30:end])[1]) |> pmean
    100perf + robust + 0.002noise
end

cost(params)

## We can now perform the optimization and visualization again

res2 = optimize(cost, params, Optim.Options(
    show_trace        = true,
    show_every        = 50,
))
plot_optimized(params, res2)
# We should now hopefully see that the constraint on Ms is fulfilled for all realizations of the uncertain system.
res.minimizer
res2.minimizer