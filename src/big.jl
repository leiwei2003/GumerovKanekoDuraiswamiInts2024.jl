using GumerovKanekoDuraiswamiInts2024

h4 = 100.0
h3 = 0.0
h2 = 0.0
#h1 = sqrt(2)/2000
h1 = 4.9065389333867974e-18
P1 = 0.01131370849898476

I0 = GumerovKanekoDuraiswamiInts2024.I0(P1, h1, h2, h3, h4)

zerotol = 3e-16

P = P1
##
P = BigFloat(P)
h1 = BigFloat(h1)
h2 = BigFloat(h2)
h3 = BigFloat(h3)
h4 = BigFloat(h4)
##
hh = h1^2 + h2^2 + h3^2 + h4^2
R = sqrt(P^2 + hh)
Phi1 = 0.5 * log((P + R)^2 / hh) / P
R1 = sqrt(P^2 + h1^2)
hh1 = hh - h1^2
zero2 = zerotol * zerotol


##
if h1 * P < zerotol
    Phi2 = 1 / (hh + h4 * R) # catch P = 0
else
    Phi2 = atan(h1 * P / (hh + h4 * R)) / (h1 * P)
end
Phi3 = 0.5 / R1 * log((R1 + R)^2 / hh1)
I = (
    (h1^2 - 3 * h4^2) * Phi1 - h4 * (3 * h1^2 - h4^2) * Phi2
    + 3 * h4^2 * Phi3 - h4^2 / (R + h4)
    ) / (6 * h1^2) # Case 6