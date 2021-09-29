using ControlSystems, CUDA
A   = cu([-1.2 0; 0.1 -0.4]);  B = cu([1.2; 0])
C   = cu([1 0]);               D = cu([0])
sys = HeteroStateSpace(A, B, C, D)
x0  = cu(zeros(2, 1))
u(x,t)  = cu(ones(1, 1))
y, t, x = lsim(sys, u, 0:0.01:5; x0)
