using GumerovKanekoDuraiswamiInts2024

h4 = 100.0
h3 = 0.0
h2 = 0.0
h1 = sqrt(2)/2000
h1 = 1e-9
P1 = 0.01

P = P1

I0 = GumerovKanekoDuraiswamiInts2024.I0(P1, h1, h2, h3, h4)

zerotol = 3e-16

#=
P = BigFloat(P)
h1 = BigFloat(h1)
h2 = BigFloat(h2)
h3 = BigFloat(h3)
h4 = BigFloat(h4)
=#

hh = h1^2 + h2^2 + h3^2 + h4^2
h = sqrt(h1^2 + h2^2 + h3^2 + h4^2)
R = sqrt(P^2 + hh)
Phi1 = 0.5 * log((P + R)^2 / hh) / P
R1 = sqrt(P^2 + h1^2)
hh1 = h2^2 + h3^2 + h4^2
zero2 = zerotol * zerotol




#if h1 * P < zerotol
#    Phi2 = 1 / (hh + h4 * R) # catch P = 0
#else
    Phi2 = atan(h1 * P / (hh + h4 * R)) / (h1 * P)
#end
Phi3 = 0.5 / R1 * log((R1 + R)^2 / hh1)
I = (
    (h1^2 - 3 * h4^2) * Phi1 - h4 * (3 * h1^2 - h4^2) * Phi2
    + 3 * h4^2 * Phi3 - h4^2 / (R + h4)
    ) / (6 * h1^2) # Case 6


I = Phi1/6 + 1/2 * (h4/h1)^2 * (Phi3 - Phi1) + 1/6*(h4/h1)^2*(h4* Phi2 - 1/(R + h4) ) - 1/2 * h4 * Phi2 

I = Phi1/6 - 1/2 * h4 * Phi2 

(h4/h1)^2 / 6 * (((3*Phi3 + h4* Phi2) - 1/(R + h4) - 3*Phi1) )

Ref = (Phi3 - Phi1)*(Phi3 + Phi1)
Diff = P^2 * (((R1+R)^2/hh1)^2 - ((P+R)^2/hh)^2 - 2 * (R1+R)^2/hh1 + 2 * ((P+R)^2/hh))
        - h1^2 * (((P+R)^2/hh)^2 - 2 * ((P+R)^2/hh) + 1)

Diff2= 1/4 * 1/((P^2+h1^2)*P^2)*Diff

final = (1/4 * 1/((P^2+h1^2)*P^2)*(- h1^2 * (((P+R)^2/hh)^2 - 2 * ((P+R)^2/hh) + 1)))/(Phi1+Phi3)/h1^2

#

Ref = Phi3^2 - Phi1^2
Alt = 1/4 * 1/((P^2+h1^2)*P^2) * (P^2 * log((R1 + R)^2 / hh1)^2 - (P^2 + h1^2)*log((P + R)^2 / hh)^2)
Alt2 = 1/4 * 1/((P^2+h1^2)*P^2) * (P^2 * (log((R1 + R)^2 / hh1)^2 - log((P + R)^2 / hh)^2) - h1^2*log((P + R)^2 / hh)^2)

Alt3 = 1/4 * 1/((P^2+h1^2)*P^2) * ( 
    P^2 * (log((R1 + R)^2 / hh1) - log((P + R)^2 / hh)) * log((R1 + R)^2 / hh1) + log((P + R)^2 / hh) * (P^2 * log((R1 + R)^2 / hh1) - P^2*log((P + R)^2 / hh) - h1^2 * log((P + R)^2 / hh)) 
    )

a = (R1+R)^2/hh1 # in Phi3
b = (P+R)^2/hh # in Phi1

New = 1/(P^2*(h1^2+P^2)) * 1/4 * (P^2 * log(a/b) * log(a*b) - h1^2 * log(b)^2) # = Ref


Ref2 = 1/h1^2 * log(a/b) * log(a*b) - log(b)^2 / P^2 # = Ref * (P^2*(h1^2+P^2)) * 4 / (P^2 * h1^2)

Ref3 = 1/h1^2 * (a/b -1) * log(a*b) - log(b)^2 / P^2   # log(1+e) ≈ e für e<<

#
Ref4 = (a/b - 1)
Ref5 = (R1+R)^2/(P+R)^2 * hh/hh1 - 1 # Do NOT cancel hh/hh1 to 1 !!

Ref6 = 1/(h4^2 * (P + R)^2) * (h4^2 * (h1^2 + h1^2 * R/P) + h1^2 * (h4^2 + 2 * (P^2 + h1^2 + R * R1)))

Ref7 = h1^2/(h4^2 * (P + R)^2) * (h4^2 * (2 + R/P) + 2 * (h1^2 + P^2 + R1 * R))

# durch h1^2 teilen, wegen Formel von I

Ref11 = 1/(h4^2 * (P + R)^2) * (h4^2 * (2 + R/P) + 2 * (h1^2 + P^2 + R1 * R)) # BigFloat = Float64 :)

#
I1 = Phi1/6 + 1/2 * (h4/h1)^2 * (Phi3 - Phi1) + 1/6*(h4/h1)^2*(h4* Phi2 - 1/(R + h4) ) - 1/2 * h4 * Phi2 # Case 6

I2 = Phi1/6 + 1/2 * (h4/h1)^2 * (Phi3^2 - Phi1^2)/(Phi1 + Phi3) + 1/6 * (h4/h1)^2 * (h4* Phi2 - 1/(R + h4) ) - 1/2 * h4 * Phi2 # Case 6

I3 = Phi1/6 + 1/2 * (h4/h1)^2 * (
    (P^2 * log(a*b) * log(a/b) - h1^2 * log(b)^2) / (4 * P^2 * (P^2 + h1^2)) # Phi3^2 - Phi1^2
)/(Phi1 + Phi3) + 1/6 * (h4/h1)^2 * (h4* Phi2 - 1/(R + h4)) - 1/2 * h4 * Phi2
# log(a/b) ersetzen
I4 = Phi1/6 + 1/2 * (h4/h1)^2 * (
    (P^2 * log(a*b) * (a/b - 1) - h1^2 * log(b)^2) / (4 * P^2 * (P^2 + h1^2)) # Phi3^2 - Phi1^2
)/(Phi1 + Phi3) + 1/6 * (h4/h1)^2 * (h4* Phi2 - 1/(R + h4)) - 1/2 * h4 * Phi2
# a/b - 1 ersetzen
I5 = Phi1/6 + 1/2 * (h4/h1)^2 * (
    (P^2 * log(a*b) * Ref7 - h1^2 * log(b)^2) / (4 * P^2 * (P^2 + h1^2)) # Phi3^2 - Phi1^2
)/(Phi1 + Phi3) + 1/6 * (h4/h1)^2 * (h4* Phi2 - 1/(R + h4)) - 1/2 * h4 * Phi2
# h1^2 kürzen
I6 = Phi1/6 + 1/2 * h4^2 * (
    (P^2 * log(a*b) * Ref11 - log(b)^2) / (4 * P^2 * (P^2 + h1^2)) # (Phi3^2 - Phi1^2) / h1^2
)/(Phi1 + Phi3) + 1/6 * (h4/h1)^2 * (h4* Phi2 - 1/(R + h4)) - 1/2 * h4 * Phi2

#

Var = h4 * Phi2 - 1/(R + h4)

Var2 = h4 / (h1^2 + h4^2 + h4 * R) - 1/(R + h4) # atan(x) ≈ x

Var25 = -h1^2 / ((hh + R * h4) * (R + h4)) # just as precise as Var3

Var3 = -h1^2 / (2 * hh * h4 + h4 * P^2 + (2*h4^2 + h1^2) * R)

# -h1^2 mit Formel kürzen
Var4 = 1 / (2 * hh * h4 + h4 * P^2 + (2*h4^2 + h1^2) * R)

#

In1 = Phi1/6 + 1/2 * h4^2 * (
    (P^2 * log(a*b) * Ref11 - log(b)^2) / (4 * P^2 * (P^2 + h1^2)) # (Phi3^2 - Phi1^2) / h1^2
)/(Phi1 + Phi3) + 1/6 * (h4/h1)^2 * Var3 - 1/2 * h4 * Phi2

In2 = Phi1/6 + 1/2 * h4^2 * (
    (P^2 * log(a*b) * Ref11 - log(b)^2) / (4 * P^2 * (P^2 + h1^2)) # (Phi3^2 - Phi1^2) / h1^2
)/(Phi1 + Phi3) - 1/6 * h4^2 * Var4 - 1/2 * h4 * Phi2


# new method for Phi3^2 - Phi1^2

cPhi = 1/2 * ((h1/h4)^2 - 1/2 * (h1/h4)^4 + 1/3 * (h1/h4)^6)
    + 1/R * (R1 - P - 1/2 * h1^2/R + 1/3 * (R1^3 - P^3)/(P^2 + h1^2 + h4^2))

I6new = Phi1/6 + 1/2 * (h4/h1)^2 * 1/(Phi1 + Phi3) * 1/(P^2 * (P^2 + h1^2)) *
        (P^2 * (log((R1 + R)/h4) + log((P + R)/h)) * cPhi - h1^2 * log((P + R)/h)^2)
    + 1/6*(h4/h1)^2*(h4* Phi2 - 1/(R + h4) ) - 1/2 * h4 * Phi2