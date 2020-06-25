using NLsolve

λ_p    = 0.15             # net markup
α      = 1-0.60/(1+λ_p)
δ    = 0.025
π_ss = 1.02^(1/4)
γ    = log(1.02)/4
disc = 0.99

h    = 0.5
g_x  = 0.18

s    = 0.995
θ    = 1/4
f_s1 = 3/4
f_h1 = 1/4
A2   = 0.7
v    = 0.7
t_x  = 0.2
τ_x  = 0
#ψ_h1 = θ * f_h1 / 1.1
#ψ_h2 = θ * (1 - f_h1) / 1.1
ψ_h1 = 0
ψ_h2 = 0
H1   = 1
H2   = 1


L1   = ((1-θ) * f_s1 + θ * f_h1) * H1
L2   = ((1-θ) * (1 - f_s1) + θ * (1 - f_h1)) * H2

ρ    = exp(γ) / disc - (1 - δ)

w1   = (α ^ α * (1 - α) ^ (1 - α) / (1 + λ_p) / (ρ ^ α)) ^ (1 / (1 - α))
w2   = A2 * w1

klr1 = w1 / ρ * α / (1 - α)
klr2 = A2 * klr1
k    = klr1 * L1

F1   = λ_p / (1 + λ_p) * L1 * klr1 ^ α
y1   = L1 * klr1 ^ α - F1
y    = y1 / v
y2   = (1 - v) * y
F2   = λ_p * y2 / A2
x    = y

is   = (1 - (1 - δ) * exp(-γ)) * exp(γ) * k / (1 - θ)
g    = g_x * x
c    = y - g - (1 - θ) * is

τ    = τ_x * x
τ_h1 = τ * ψ_h1 / θ / (1 - f_h1)
τ_h2 = τ * ψ_h2 / θ / (1 - f_h1)
τ_s  = τ - θ * (f_h1 * τ_h1 + (1 - f_h1) * τ_h2)

t    = t_x * x
t_h1 = t * ψ_h1 / θ / f_h1
t_h2 = t * ψ_h2 / θ / (1 - f_h1)
t_s  = (t - θ * ( f_h1 * t_h1 + (1 - f_h1) * t_h2)) / (1 - θ)

c_h10 = c / θ / f_h1
c_h20 = c / θ / (1 - f_h1)
c_s0  = c / (1 - θ)
λ_h10 = exp(γ) / (exp(γ) * c_h10 - h * c)
λ_h20 = exp(γ) / (exp(γ) * c_h20 - h * c)
λ_s0  = exp(γ) / (exp(γ) * c_s0 - h * c)
bR0   = (g - t + τ) / (1 - (1.025 / 4) / exp(γ) / π_ss)
R0    = (1.025 / 4)

function distSS(x)
global w1, w2, H1, H2, th1, th2, τ_h1, τ_h2, c, g, t, τ
global s, θ, f_h1, f_s1, γ, π_ss, h, disc

c_h1 = x[1]
c_h2 = x[2]
c_s   = x[3]
λ_h1 = x[4]
λ_h2 = x[5]
λ_s  = x[6]
bR   = x[7]
R    = x[8]

r = zeros(8)

r[1] = w1 * H1 - t_h1 + τ_h1 + R * bR * f_s1 * (1 - s) / θ / f_h1 / exp(γ) / π_ss - c_h1
r[2] = w2 * H2 - t_h2 + τ_h2 + R * bR * (1-f_s1) * (1 - s) / θ / (1 - f_h1) / exp(γ) / π_ss - c_h2
r[3] = c / (1 - θ) - θ * (f_h1 * c_h1 + (1 - f_h1) * c_h2) / (1 - θ) - c_s

r[4] = exp(γ) / (exp(γ) - h) / c_h1 - λ_h1
r[5] = exp(γ) / (exp(γ) - h) / c_h2 - λ_h2
r[6] = exp(γ) / (exp(γ) - h) / c_s  - λ_s

r[7] = (g - t + τ) / (1 - R / exp(γ) / π_ss) - bR

r[8] = exp(γ) * π_ss / disc * λ_s / (s * λ_s + (1 - s) * (f_s1 * λ_h1 + (1 - f_s1) * λ_h2)) - R

return r

end


init = [c_h10 c_h20 c_s0 λ_h10 λ_h20 λ_s0 bR0 R0]
ss = nlsolve(distSS, init)

solution = ss.zero

println(solution)
