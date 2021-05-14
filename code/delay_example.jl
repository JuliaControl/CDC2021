# Analysis of Smith predictor + PI controller for controlling a first-order system with a delay.
# This script replicates some of the results of Ch. 8 in "Advanced PID control" by Åström & Hägglund (2006)

using ControlSystems
using Plots
using Plots.Measures

𝕗 = (14, "Computer Modern")
gr(titlefont=𝕗, tickfont=𝕗, legendfont=𝕗, guidefont=𝕗, grid=false, legend=false, bottommargin=2mm)

## Set up the nominal system and PI controller
P0 = ss(-1, 1, 1, 0)

# PI controller for nominal system P0
# To verify the pole placement, use, e.g., dampreport(feedback(P0, C0))
ω0 = 2
ζ = 0.7
C0 = placePI(P0, ω0, ζ)

## Setup delayed plant + Smith predictor-based
#  controller for a given delay τ
τ = 8
P = delay(τ) * P0
C_sp = feedback(C0, (1.0 - delay(τ))*P0)

## Plot the closed loop response 
# Reference step at t = 0 and load disturbance step at t = 15
G = [feedback(P*C_sp, 1) feedback(P, C_sp)]
fig_timeresp = lsimplot(G, t -> [1; t >= 15], 0:0.1:40,  title="τ = $τ", size=(400,270))

## Plot the frequency response of the predictor part and compare to a negative delay
C_pred = feedback(1, C0*(ss(1.0) - delay(τ))*P0)
fig_bode = bodeplot([C_pred, delay(-τ)], exp10.(-1:0.002:0.4), ls=[:solid :solid :dash :dash], title="")
plot!(yticks=[0.1, 1, 10], sp=1)
plot!(yticks=0:180:1080, sp=2)

## Check the Nyquist plot
# Note that the Nyquist curve encircles -1 for τ > 2.99
fig_nyquist = nyquistplot(C_sp * P, exp10.(-1:1e-4:2), size=(400,270), title="τ = $τ")


## Save figures
savefig(fig_timeresp, "ex1_timeresp_tau$τ.pdf")
savefig(fig_nyquist, "ex1_nyquist_tau$τ.pdf")
savefig(fig_bode, "ex1_bode_tau$τ.pdf")
