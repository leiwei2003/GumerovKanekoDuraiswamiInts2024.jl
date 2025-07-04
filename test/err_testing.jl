using GumerovKanekoDuraiswamiInts2024

function I6jl(P, h1, h4)

    hh = h1^2 + h4^2
    hh1 = h4^2
    R = sqrt(P^2 + hh)
    R1 = sqrt(P^2 + h1^2)
    
    Phi1 = 1/2 * log((P + R)^2 / hh) / P
    Phi2 = 1 / (hh + h4 * R)
    Phi3 = 1/2 / R1 * log((R1 + R)^2 / hh1)

    I = Phi1/6 + 1/2 * h4^2 * (
        (P^2 * log((P + R)^2 * (R1 + R)^2 / (hh * hh1)) * 
        1 / (h4^2 * (P + R)^2) * (h4^2 * (2 + R/P) + 2 * (h1^2 + P^2 + R1 * R)) - log((P + R)^2 / hh)^2) / (4 * P^2 * (P^2 + h1^2))
    )/(Phi1 + Phi3) - 1/6 * h4^2 * 1 / (2 * hh * h4 + h4 * P^2 + (2*h4^2 + h1^2) * R) - 1/2 * h4 * Phi2

    return I
end

h4 = 100.0
h3 = 0.0
h2 = 0.0
#h1 = sqrt(2)/2000
h1 = 4.9065389333867974e-18
P1 = 0.01131370849898476

I = I6jl(Float64(P), Float64(h1), Float64(h4))

for i in 1:50
    h1 = 2^i*1e-18
    h4 = 100.0
    P = P1

    h1 = BigFloat(h1)
    h4 = BigFloat(h4)
    P = BigFloat(P)
 

    println("\nP = ", P, ", h4 = ", h4, ", h1 = ", h1)
    println("I = ", I6jl(P, h1, h4))
end