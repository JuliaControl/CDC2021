using MonteCarloMeasurements, Optim, ControlSystems
±(m,s) = m + s*Particles(200)
unsafe_comparisons(true)
p = 1 ± 0.1
ζ = 0.3 ± 0.05
ω = 1 ± 0.05
const P = tf([p*ω], [1, 2ζ*ω, ω^2]) |> ss
const Ω = exp10.(LinRange(-2,3,100))
params = log.([1,0.1,0.1]) # Initial guess
const Msc = 1.2 # Constraint on Ms

function systems(params)
    kp,ki,kd = exp.(params)
    C  = ss(pid(; kp, ki, kd)*tf(1, [0.05, 1])^2)
    G  = feedback(P*C) # Closed-loop system
    S  = 1/(1 + P*C)   # Sensitivity function
    CS = C*S           # Noise amplification
    Gd = c2d(G, 0.1)   # Discretize the system
    y,t = vec.(step(Gd,15)) # Step-response
    C, G, S, CS, y, t
end

function cost(params::AbstractVector{T}) where T
    C, G, S, CS, y, t = systems(params)
    # Maximum of the sensitivity function
    Ms = maximum(bode(S, Ω)[1]) 
    q  = quantile(Ms, 0.9)
    # This is our performance measure
    performance = mean(abs, 1 .- y)   
    # This is our robustness constraint
    robustness = (q > Msc ? 10000(q-Msc) : zero(T)) 
    # Penalty for high variance in performance
    variance = std(performance)     
    noise = mean(sum(bode(CS, Ω[end-30:end])[1]))
    100mean(performance) + robustness +
              10variance + 0.002noise
end

res = optimize(cost, params)

## We can now perform the same computations as above to visualize the found controller
using Plots.Measures
𝕗 = (14, "Computer Modern")
gr(titlefont=𝕗, tickfont=𝕗, legendfont=𝕗, guidefont=𝕗, grid=false)
fig = plot(layout=2, size=(1000,400), bottommargin=2mm)
for (i,params) = enumerate((params, res.minimizer))
    C, G, S, CS, y, t = systems(params)
    mag   = bode(S, Ω)[1][:]
    plot!(t,y[:], title="Time response", subplot=1, lab= i==1 ? "Initial" : "Optimized", legend=:bottomright)
    plot!(Ω, mag, title="Sensitivity function", xscale=:log10, yscale=:log10, subplot=2, lab= i==1 ? "Initial" : "Optimized", legend=:bottomright)
end
hline!([Msc], l=(:black, :dash), subplot=2, lab="")
display(fig)
