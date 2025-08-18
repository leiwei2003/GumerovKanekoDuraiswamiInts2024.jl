
h2 = 1e-9
h4 = 100
P = 0.01
#=
h2 = BigFloat(h2)
h4 = BigFloat(h4)
P = BigFloat(P)
=#
h = sqrt(h2^2 + h4^2)
hh = h2^2 + h4^2
R = sqrt(P^2 + hh)
R2 = sqrt(P^2 + h2^2)
R4 = sqrt(P^2 + h4^2)
RR = P^2 + hh
RR2 = P^2 + h2^2
RR4 = P^2 + h4^2
    
Phi1 = 1/2 / P * log((P + R)^2 / hh)
Phi2 = atan(h2 * P / (hh + h4 * R)) / P
Phi4 = h4^2 / (h2 * P^2) * ((R2 / h2 * log((R2 + R) / h4) - log((h2 + h) / h4)))

I1 = (
    (1 + 3 * (h4/h2)^2) * Phi1 - 2 * (h4/h2)^3 * Phi2 - 3 * Phi4 +
    (2 * h4^2 - h2^2)/(h2^2 * (R+h))
    )/6

# nur L'Hospital (ist falsch, felht die Hälfte. mache grad neu)

I2 = (
    6/(R+h) - (6*h4^2 + h2^2)/(RR*h + R*hh)^2 * (2/(RR*h + R*hh) + (2*h4^2 - h2^2)*(2*h + 2*RR + 2*R + hh/R))
    - 2 * h4^3/((hh+h4*R)^2 + h2^2*P^2) * (
        - 2 * h2^2 * (hh + h4*R)*(2 + h4/R)
        + (hh + h4*R)^2 * ((2 + h4/R) - h2*(2 - 2*h2*h4/R))
        - 2 * ((hh + h4*R)*(2+h4/R) + P^2) * ((hh+h4*R)^2 - h2^2*(2+h4/R))
    )
)/6

# auch h2 -> 0

I3 = 1/(R4+h4) - h4^2/(RR4*h4 + R4*h4^2)^2 * (2/(RR4*h4 + R4*h4^2) + 2*h4^2 * (2*h4 + 2*RR4 + 2*R4 + h4^2/R4))
    - 1/3 * h4^3 * ((2+h4/R4) - 2 * ((h4^2 + h4*R4)*(2+h4*R4) + P^2))