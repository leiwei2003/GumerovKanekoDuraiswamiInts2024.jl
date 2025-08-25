
h1 = 80
h4 = 100
P = 0.008
#=
h1 = BigFloat(h1)
h4 = BigFloat(h4)
P = BigFloat(P)
=#
hh = h1^2 + h4^2
hh1 = h4^2
R = sqrt(P^2 + hh)
R1 = sqrt(P^2 + h1^2)
R4 = sqrt(P^2 + h4^2)
RR = P^2 + hh
RR1 = P^2 + h1^2
RR4 = P^2 + h4^2
    
Phi1 = 1/2 / P * log((P + R)^2 / hh)
Phi2 = atan(h1 * P / (hh + h4 * R)) / (h1 * P)
Phi3 = 1/2 / R1 * log((R1 + R)^2 / hh1)

I = Phi1/6 + 1/2 * (h4/h1)^2 * (Phi3 - Phi1) + 1/6 * (h4/h1)^2 * (h4 * Phi2 - 1/(R + h4)) - 1/2 * h4 * Phi2

# 1/h1^2 * (Phi3 - Phi1)

Ln1 = 1/h1^2 * (Phi3 - Phi1)

Ln2 = 1/(2 * P * R1) * (P / (R1 * R) - 1/2 * 1/R1 * log((P + R)^2/hh) - R1 * (1/(R * (P + R)) - 1/hh)) # bis h1 < 1e-6

Ln3 = 1/(2 * P^2) * (1/R4 - 1/P * log((P + R4)/h4) - P * (1/(RR4 + P * R4) - 1/h4^2)) # bis h1 < 1e-9

# 1/h1^2 * (h4 * Phi2 - 1/(R + h4))

Tan1 = 1/h1^2 * (h4 * Phi2 - 1/(R + h4))

Tan2 = -h4/3 * 1/((h1 * P)^2 + (hh + h4 * R)^2)^2 * (
    (hh + h4 * R)^2 + P^2 * (h4^2 + h4 * R) + 2 * ((h4^2 + h4 * R)^2 - h1^4)) # bis h1 < 1e-3, sehr gut

Tan3 = -h4/3 * 1/(h4^2 + h4 * R4)^3 * (3 * (h4^2 + h4 * R4) + P) # bis h1 < 1e-3, konst. Fehler

# I6

I1 = Phi1/6 + 1/2 * (h4/h1)^2 * (Phi3 - Phi1) + 1/6 * (h4/h1)^2 * (h4 * Phi2 - 1/(R + h4)) - 1/2 * h4 * Phi2

# nicht weiter gekürzt, nur L'Hopital
I2 = Phi1/6 + 1/4 * h4^2/RR1 * (
        1/R - Phi1 - RR1/P * (1/(R * P + RR) - 1/hh)
    ) - 1/18 * h4^3/((h1 * P)^2 + (hh + h4 * R)^2)^2 * (
        (hh + h4 * R)^2 + P^2 * (h4^2 + h4 * R) + 2 * ((h4^2 + h4 * R)^2 - h1^4)
    ) - 1/2 * h4 * Phi2

# alles gekürzt
I3 = 1/6 * log((P + R4)/h4) / P + 1/4 * (h4/P)^2 * (
        1/R4 - 1/P * (log((P + R4)/h4)) - P * (1/(RR4 + P * R4) - 1/h4^2)
    ) - 1/18 * 1/(h4 + R4)^3 * (3 * (h4^2 + h4 * R4) + P) - 1/2 * 1/(h4 + R4)
