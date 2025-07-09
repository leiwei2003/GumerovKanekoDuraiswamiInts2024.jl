
h4 = 100
h1 = sqrt(2)/2000
P = 0.01



function taylorSqrt(x, n)

    sum = 0
    for i in 0:n
        sum += 2 / (i + 1) * binomial(2 * i, i) * (-x/4)^(i+1)
    end

    return 1 - sum
end

function I_APRX(P, h1, h4, n)

    hh = h1^2 + h4^2
    R = sqrt(P^2 + hh)
    R1 = sqrt(P^2 + h1^2)
    hh1 = h4^2

    Phi1 = 0.5 * log((P + R)^2 / hh) / P
    Phi2 = atan(h1 * P / (hh + h4 * R)) / (h1 * P)
    Phi3 = 0.5 / R1 * log((R1 + R)^2 / hh1)

    a = (R1+R)^2/hh1 # in Phi3
    b = (P+R)^2/hh # in Phi1

    Tayl = h4^2 * (h1^2 + 2 * P * R * (taylorSqrt((h1/P)^2, n)) - 1) + h1^2 * (2 * P^2 + 2 * h1^2 + h4^2 + 2 * R * R1)

    Ref11 = 1/(h4^2 * (P + R)^2) * Tayl
    Var4 = 1 / (2 * hh * h4 + h4 * P^2 + (2*h4^2 + h1^2) * R)

    return Phi1/6 + 1/2 * h4^2 * (
                (P^2 * log(a*b) * Ref11 - log(b)^2) / (4 * P^2 * (P^2 + h1^2)) # (Phi3^2 - Phi1^2) / h1^2
            )/(Phi1 + Phi3) - 1/6 * h4^2 * Var4 - 1/2 * h4 * Phi2
end

# (h4^2 * (h1^2 + 2 * P * R * taylorSqrt((h1/P)^2, n)) + h1^2 * (2 * P^2 + 2 * h1^2 + h4^2 + 2 * R * R1))

function I_REAL(P, h1, h4)

    hh = h1^2 + h2^2 + h3^2 + h4^2
    R = sqrt(P^2 + hh)
    R1 = sqrt(P^2 + h1^2)
    hh1 = h2^2 + h3^2 + h4^2

    Phi1 = 0.5 * log((P + R)^2 / hh) / P
    Phi2 = atan(h1 * P / (hh + h4 * R)) / (h1 * P)
    Phi3 = 0.5 / R1 * log((R1 + R)^2 / hh1)

    return ((h1^2 - 3 * h4^2) * Phi1 - h4 * (3 * h1^2 - h4^2) * Phi2+ 3 * h4^2 * Phi3 - h4^2 / (R + h4)) / (6 * h1^2)
end

h1 = BigFloat(h1)
h4 = BigFloat(h4)
P = BigFloat(P)

for k in 0:10
    println("Order = ", k + 1, ", Res = ", taylorSqrt((h1/P)^2, k))
end

for k in 0:10
    println("Order = ", k + 1, ", Err = ", abs(taylorSqrt((h1/P)^2, k)-sqrt(1+(h1/P)^2))/sqrt(1+(h1/P)^2))
end

