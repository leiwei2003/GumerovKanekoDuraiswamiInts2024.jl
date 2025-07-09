using GumerovKanekoDuraiswamiInts2024

h4 = 100.0
h3 = 0.0
h2 = 0.0
h1 = sqrt(2)/2000
#h1 = 4.9065389333867974e-18
P1 = 0.01

P = P1

I0 = GumerovKanekoDuraiswamiInts2024.I0(P1, h1, h2, h3, h4)

zerotol = 3e-16

#
P = BigFloat(P)
h1 = BigFloat(h1)
h2 = BigFloat(h2)
h3 = BigFloat(h3)
h4 = BigFloat(h4)
#
I1 = GumerovKanekoDuraiswamiInts2024.I0(P, h1, h2, h3, h4)


hh = h1^2 + h2^2 + h3^2 + h4^2
R = sqrt(P^2 + hh)
Phi1 = 0.5 * log((P + R)^2 / hh) / P
R1 = sqrt(P^2 + h1^2)
hh1 = h2^2 + h3^2 + h4^2
zero2 = zerotol * zerotol

Phi1 = 0.5 * log((P + R)^2 / hh) / P
Phi2 = atan(h1 * P / (hh + h4 * R)) / (h1 * P)
Phi3 = 0.5 / R1 * log((R1 + R)^2 / hh1)

I = (
    (h1^2 - 3 * h4^2) * Phi1 - h4 * (3 * h1^2 - h4^2) * Phi2
    + 3 * h4^2 * Phi3 - h4^2 / (R + h4)
    ) / (6 * h1^2) # Case 6

I2 = Phi1/6 + 1/2 * (h4/h1)^2 * (Phi3 - Phi1) + 1/6 * (h4/h1)^2 * (h4 * Phi2 - 1/(R + h4)) - 1/2 * h4 * Phi2 

I6 = Phi1/6 + 1/2 * h4^2 * (
        (P^2 * log((P + R)^2 * (R1 + R)^2 / (hh * hh1)) * 
        1 / (h4^2 * (P + R)^2) * (h4^2 * (2 + R/P) + 2 * (h1^2 + P^2 + R1 * R)) - log((P + R)^2 / hh)^2 ) / (4 * P^2 * (P^2 + h1^2))
        )/(Phi1 + Phi3) - 1/6 * h4^2 * 1 / (2 * hh * h4 + h4 * P^2 + (2*h4^2 + h1^2) * R) - 1/2 * h4 * Phi2


test = 3*h4^2/h1^2 * (Phi3 - Phi1)

test2 = 3*h4^2/(2*P*h1^2) * log((R1+R)^2/(P+R)^2 * (h1^2+h4^2)/h4^2)

test3 = 3*h4^2/(2*P*h1^2) * ((R1+R)^2/(P+R)^2 * (h1^2+h4^2)/h4^2 - 1)




