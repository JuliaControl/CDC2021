using ControlSystems, Plots

P₀ = tf(1.0, [1, 1])
ω₀ = 2; ζ = 0.7
kp, ki, C₀ = placePI(P₀, ω₀, ζ)
τ = 8.0
P = delay(τ) * P₀
Cₛₚ = feedback(C₀, (1 - delay(τ))*P₀)
t = 0:0.1:40
y = step_ref_disturbance(P, Cₛₚ, 0:0.1:40, 15)
plot(t, y)
C_pred = feedback(1, C₀*(ss(1.0) - delay(τ))*P₀)
bodeplot([C_pred,delay(-τ)],exp10.(-1:0.002:1))	
nyquistplot(Cₛₚ*P, exp10.(-1:0.0001:1))