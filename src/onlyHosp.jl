h1 = 0
h2 = 1e-3
h3 = 0
h4 = 100
P = 0.01
#=
h1 = BigFloat(h1)
h2 = BigFloat(h2)
h3 = BigFloat(h3)
h4 = BigFloat(h4)
P = BigFloat(P)
=#
h = sqrt(h1^2 + h2^2 + h3^2 + h4^2)
hh = h1^2 + h2^2 + h3^2 + h4^2
R = sqrt(P^2 + hh)
R1 = sqrt(P^2 + h1^2)
R2 = sqrt(P^2 + h2^2)
R4 = sqrt(P^2 + h4^2)
RR = P^2 + hh
RR1 = P^2 + h1^2
RR2 = P^2 + h2^2
RR4 = P^2 + h4^2

# erster Teil
K1 = 3/P * (h4/h2)^2 * log((P + R)/h) - 3/P^2 * (h4/h2)^2 * R2 * log((R2 + R)/h4) + 3/P^2 * h4^2/h2 * log((h2 + h)/h4)
K2 = - 3/2 * (h4/P)^2 * (1/(P + R) + P/hh + 1/R2 * log((R2 + R)/h4) - 1/h - 1/h2 * log((h2 + h)/h4))
K3 = - 3/2 * (h4/P)^2 * (1/(P + R) + P/hh + 1/R2 * log((R2 + R)/h4) - 2/h)