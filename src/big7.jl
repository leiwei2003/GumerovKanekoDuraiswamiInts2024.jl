
h2 = 1000000
h4 = 100
P = 0.008
#
h2 = BigFloat(h2)
h4 = BigFloat(h4)
P = BigFloat(P)
#
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
    (2 * h4^2 - h2^2)/(h2^2*(R + h))
    )/6

# h2 -> 0

I2 = (
    1/P * log((P + R4)/h4) # Phi1 auch gut :)
    - 3/2 * h4^2/P^2 * (1/(P + R4) + P/h4^2 + 1/P * log((P + R4)/h4) - 2/h4) # yay
    - 1/(R4 + h4) * (1 + h4/3 * (1/R4 - 1/(R4 + h4) * (4 + h4/R4 + 2/(R4 + h4) * P^2/h4))) # Phi2 gut :)
    )/6

# nur L' Hospital

I3 = 1/6 * (
    1/P * log((P + R)/h)
    - 3/2 * (h4/P)^2 * (1/(P + R) + P/hh + 1/R2 * log((R2 + R)/h4) - 2/h)
    + 1/6 * ((6*h4^2 - 2*h2^2 - R*h)/(RR*h + R*hh)
        - 4*h4^3*(hh + h4*R)/((hh + h4*R)^2 + (h2*P)^2)^2 * (P^2 + (2 + h4/R)*(hh + h4*R)))
    )

# Phi1 u. Phi4
Test1 = 3 *((h4/h2)^2 * Phi1 - Phi4)

# L' Hospital
Test2 = 3/2 * h4^2/P^2 * (2/h - 1/(P + R) - 1/R2 * log((R2 + R)/h4))