using GumerovKanekoDuraiswamiInts2024

function I6jl(P, h1, h4)

    hh = h1^2 + h4^2
    hh1 = h4^2
    R = sqrt(P^2 + hh)
    R1 = sqrt(P^2 + h1^2)
    
    Phi1 = 1/2 * log((P + R)^2 / hh) / P
    Phi3 = 1/2 / R1 * log((R1 + R)^2 / hh1)

    if true#h1 < 1e-6
        Phi2 = 1 / (hh + h4 * R)

        I = Phi1/6 + 1/2 * h4^2 * (
            (P^2 * log((P + R)^2 * (R1 + R)^2 / (hh * hh1)) * 
            1 / (h4^2 * (P + R)^2) * (h4^2 * (2 + R/P) + 2 * (h1^2 + P^2 + R1 * R)) - log((P + R)^2 / hh)^2) / (4 * P^2 * (P^2 + h1^2))
        )/(Phi1 + Phi3) - 1/6 * h4^2 * 1 / (2 * hh * h4 + h4 * P^2 + (2*h4^2 + h1^2) * R) - 1/2 * h4 * Phi2
    else
        Phi2 = atan(h1 * P / (hh + h4 * R)) / (h1 * P)
        
        I = (
            (h1^2 - 3 * h4^2) * Phi1 - h4 * (3 * h1^2 - h4^2) * Phi2
            + 3 * h4^2 * Phi3 - h4^2 / (R + h4)
        ) / (6 * h1^2)
    end

    return I
end

function I6jlNewer(P, h1, h4)

    hh = h1^2 + h4^2
    hh1 = h4^2
    R = sqrt(P^2 + hh)
    R1 = sqrt(P^2 + h1^2)
    RR = P^2 + hh
    RR1 = P^2 + h1^2
    
    Phi1 = 1/2 / P * log((P + R)^2 / hh)
    Phi2 = atan(h1 * P / (hh + h4 * R)) / (h1 * P)
    Phi3 = 1/2 / R1 * log((R1 + R)^2 / hh1)

    Phi1/6 + 1/4 * h4^2/RR1 * (
        1/R - Phi1 - RR1/P * (1/(R * P + RR) - 1/hh)
    ) - 1/18 * h4^3/((h1 * P)^2 + (hh + h4 * R)^2)^2 * (
        (hh + h4 * R)^2 + P^2 * (h4^2 + h4 * R) + 2 * ((h4^2 + h4 * R)^2 - h1^4)
    ) - 1/2 * h4 * Phi2

    return I
end

function I6jlNewest(P, h1, h4)

    hh = h1^2 + h4^2
    hh1 = h4^2
    R = sqrt(P^2 + hh)
    R1 = sqrt(P^2 + h1^2)
    R4 = sqrt(P^2 + h4^2)
    RR = P^2 + hh
    RR1 = P^2 + h1^2
    RR4 = P^2 + h4^2
    
    Phi1 = 1/2 * log((P + R)^2 / hh) / P
    Phi3 = 1/2 * log((R1 + R)^2 / hh1) / R1
    Phi2 = 1 / (hh + h4 * R)

    I = 1/6 * log((P + R4)/h4) / P + 1/4 * (h4/P)^2 * (
            1/R4 - 1/P * (log((P + R4)/h4)) - P * (1/(RR4 + P * R4) - 1/h4^2)
        ) - 1/18 * 1/(h4 + R4)^3 * (3 * (h4^2 + h4 * R4) + P) - 1/2 * 1/(h4 + R4)

    return I
end

function I6jlOld(P, h1, h4)

    hh = h1^2 + h4^2
    hh1 = h4^2
    R = sqrt(P^2 + hh)
    R1 = sqrt(P^2 + h1^2)
    
    Phi1 = 1/2 * log((P + R)^2 / hh) / P
    Phi3 = 1/2 / R1 * log((R1 + R)^2 / hh1)

    if h1 * P < 3e-16
        Phi2 = 1 / (hh + h4 * R) # catch P = 0
    else
        Phi2 = atan(h1 * P / (hh + h4 * R)) / (h1 * P)
    end
    I = (
        (h1^2 - 3 * h4^2) * Phi1 - h4 * (3 * h1^2 - h4^2) * Phi2
        + 3 * h4^2 * Phi3 - h4^2 / (R + h4)
        ) / (6 * h1^2)

    return I
end

h4 = 100.0
h3 = 0.0
h2 = sqrt(2)/2000
h1 = 0.0
P = 0.8

I = I6jl(Float64(P), Float64(h1), Float64(h4))

function I7jlOld(h2, h4, P)

    hh = h2^2 + h4^2
    h = sqrt(h2^2 + h4^2)
    R = sqrt(P^2 + hh)
    R2 = sqrt(P^2 + h2^2)
    
    Phi1 = 1/2 * log((P + R)^2 / hh) / P
    Phi4 = h4^2 / (h2 * P^2) * ((R2 / h2 * log((R2 + R) / h4) - log((h2 + h) / h4)))
    Phi2 = atan(h2 * P / (hh + h4 * R)) / P

    return (
            (1 + 3 * h4^2 / h2^2) * Phi1 - 2 * (h4 / h2)^3 * Phi2 - 3 * Phi4 +
            (2 * h4^2 - h2^2) / (h2^2 * (R + h))
            ) / 6
end

function I7jlNew(h4, P)

    R4 = sqrt(P^2 + h4^2)

    return (
            1/P * log((P + R4)/h4) # Phi1 auch gut :)
            - 3/2 * h4^2/P^2 * (1/(P + R4) + P/h4^2 + 1/P * log((P + R4)/h4) - 2/h4) # yay
            - 1/(R4 + h4) * (1 + h4/3 * (1/R4 - 1/(R4 + h4) * (4 + h4/R4 + 2/(R4 + h4) * P^2/h4))) # Phi2 gut :)
            )/6
end

println("+++ Begin +++\n")

for i in 0:25
    h2 = 2^i*1e-5
    h4 = h4
    P = P
    #
    h2 = BigFloat(h2)
    h4 = BigFloat(h4)
    P = BigFloat(P)
    #

    #println("\nP = ", P, ", h4 = ", h4, ", h1 = ", h1, " I = ", I7jlOld(P, h2, h4))
    println("\nh2 = ", h2, "\nI = ", I7jlOld(h2, h4, P))

end


println("\n\n", I7jlOld(sqrt(2)/2000, h4, P))

#GumerovKanekoDuraiswamiInts2024.I0(0.01, sqrt(2)/2000, 0, 0, h4)