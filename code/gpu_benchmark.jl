using ControlSystems, CUDA, Plots, BenchmarkTools

function bench_sys_size(nx)
    nu = 10
    ny = 2
    dt = 0.01
    A  = randn(Float32, nx,nx)
    B  = randn(Float32, nx,nu)
    C = randn(Float32, ny,nx)
    D = zeros(Float32, ny,nu)
    sys = ss(A, B, C, D)
    sysd = c2d(sys, dt)
    sysg = HeteroStateSpace(cu(A), cu(B), cu(C), cu(D))
    sysdg = HeteroStateSpace(cu(sysd.A), cu(sysd.B), cu(sysd.C), cu(sysd.D), dt)

    tmp = ones(nu)
    u(x, t) = tmp
    tmpg = cu(ones(nu))
    ug(x, t) = tmpg

    x0 = zeros(nx)
    x0g = cu(zeros(nx))

    t = 0:dt:1

    uu = reduce(hcat, u(x0, tt) for tt in t) |> Array
    uug = reduce(hcat, ug(x0g, tg) for tg in t) |> cu

    t1 = @benchmark lsim($(sys), $(u), $(t); x0=$(x0))
    t2 = @benchmark CUDA.@sync lsim($(sysg), $(ug), $(t); x0=$(x0g))
    t3 = @benchmark lsim($(sysd), $(uu), $(t); x0=$(x0))
    t4 = @benchmark CUDA.@sync lsim($(sysdg), $(uug), $(t); x0=$(x0g))

    # Find minimum time for each test and convert to milliseconds
    return map(x->time(minimum(x)) ./ 1e6, [t1, t2, t3, t4])
end

sizes = 10:10:1000
durs = Matrix{Float64}(undef, length(sizes), 4)
for (i, n) in enumerate(sizes)
    durs[i, :] .= bench_sys_size(n)
end

𝕗 = (14, "Computer Modern")
gr(titlefont=𝕗, tickfont=𝕗, legendfont=𝕗, guidefont=𝕗, grid=false)
fig = plot(layout=(2, 1), size=(700, 500))
plot!(fig, sizes, durs[:, 3:4], title="Discrete system with array input", subplot=1, legend=:topleft, labels=["CPU" "GPU"], ylabel="Time [ms]")
plot!(fig, sizes, durs[:, 1:2], title="Continuous system with function input", subplot=2, legend=false, xlabel="State space size", ylabel="Time [ms]")
savefig("lsim_speed.pdf")
