using MonteCarloMeasurements, Optim, ControlSystems
±(m,s) = m + s*Particles(200)
unsafe_comparisons(true)
p = 1 ± 0.1
ζ = 0.3 ± 0.05
ω = 1 ± 0.05
const P = tf([p*ω], [1, 2ζ*ω, ω^2]) |> ss
const w = exp10.(LinRange(-2,3,100))
params = log.([1,0.1,0.1]) # Initial guess
const Msc = 1.2 # Constraint on Ms

function systems(params::AbstractVector{T}) where T
    kp,ki,kd = exp.(params)
    C  = ss(pid(; kp, ki, kd)*tf(1, [0.05, 1])^2)
    G  = feedback(P*C) # Closed-loop system
    S  = 1/(1 + P*C)   # Sensitivity function
    CS = C*S           # Noise amplification
    Gd = c2d(G,0.1)    # Discretize the system
    y,t = vec.(step(Gd,15)) # Step-response
    C, G, S, CS, y, t
end

function cost(params::AbstractVector{T}) where T
    C, G, S, CS, y, t = systems(params)
    # Maximum of the sensitivity function
    Ms = maximum(bode(S, w)[1]) 
    q  = quantile(Ms, 0.9)
    # This is our performance measure
    performance = mean(abs, 1 .- y)   
    # This is our robustness constraint
    robustness = (q > Msc ? 10000(q-Msc) : zero(T)) 
    # Penalty for high variance in performance
    variance = std(performance)     
    noise = mean(sum(bode(CS, w[end-30:end])[1]))
    100mean(performance) + robustness +
              10variance + 0.002noise
end

res = optimize(cost, params)