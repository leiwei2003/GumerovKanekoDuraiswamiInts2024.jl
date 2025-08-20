

function real(h1, h4, P)

    hh = h1^2 + h4^2
    hh1 = h4^2
    R = sqrt(P^2 + hh)
    R1 = sqrt(P^2 + h1^2)

    Phi1 = 1/2 / P * log((P + R)^2 / hh)
    Phi2 = atan(h1 * P / (hh + h4 * R)) / (h1 * P)
    Phi3 = 1/2 / R1 * log((R1 + R)^2 / hh1)

    return Phi1/6 + 1/2 * (h4/h1)^2 * (Phi3 - Phi1) + 1/6 * (h4/h1)^2 * (h4 * Phi2 - 1/(R + h4)) - 1/2 * h4 * Phi2
end

function aprx(h4, P)

    R4 = sqrt(P^2 + h4^2)
    RR4 = P^2 + h4^2

    return 1/6 * log((P + R4)/h4) / P + 1/4 * (h4/P)^2 * (
            1/R4 - 1/P * log((P + R4)/h4) - P * (1/(RR4 + P * R4) - 1/h4^2)
            ) - 1/18 * 1/(h4 + R4)^3 * (3 * (h4^2 + h4 * R4) + P) - 1/2 * 1/(h4 + R4)
end

H1 = BigFloat(1e-7)
H4 = BigFloat(0.001)
Pp = BigFloat(0.0001)

for h in 0:7
    for j in 0:5
        for p in 0:5

            h1 = H1 * 10^h
            h4 = H4 * 10^j
            P = Pp * 10^p

            err = abs(real(h1, h4, P) - aprx(h4, P)) / real(h1, h4, P)

            println("h1 = ", h1, " h4 = ", h4, " P = ", P, " Err = ", err)
        end
    end
end

