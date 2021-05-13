using ControlSystems, SymPy, Polynomials
s = tf("s")
@syms δ k2 γ

H = (s - 1//2) / s / (s - 1)
Gc = -k2 * (s + δ) / (s + γ)
sys_cl = feedback(H, Gc)
# Get the denominator polynomial
Q = denpoly(sys_cl)[]
# Scale so that leading coeff. equals one
Q_scaled = Q /  Q[end]
Q_desired = Polynomial([1,1])^3
SymPy.solve(coeffs(Q_scaled - Q_desired),
    (δ, k2, γ))
# Output: (-1/9, -18, -14)
