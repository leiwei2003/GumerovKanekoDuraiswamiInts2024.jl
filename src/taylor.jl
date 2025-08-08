
h4 = 1
h1 = sqrt(2)/2000
#h1 = 5e-15
P = 0.01

function taylorSqrt(x, n)

    sum = 0
    for i in 1:n
        sum += 2 / i * binomial(2 * (i - 1), i - 1) * (-x/4)^i
    end

    return 1 - sum
end

function taylorLn(x, n) # ln(x+1)

    sum = 0
    for i in 1:n
        sum -= (-1)^i * (x)^i / i
    end

    return sum
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

    Tayl = h4^2 * (h1^2 + 2 * P * R * (taylorSqrt((h1/P)^2, n) - 1)) + h1^2 * (2 * P^2 + 2 * h1^2 + h4^2 + 2 * R * R1)

    Ref6 = 1/(h4^2 * (P + R)^2) * Tayl
    Var4 = 1 / (2 * hh * h4 + h4 * P^2 + (2*h4^2 + h1^2) * R)

    return Phi1/6 + 1/2 * (h4/h1)^2 * (
                (P^2 * log(a*b) * Ref6 - h1^2 * log(b)^2) / (4 * P^2 * (P^2 + h1^2)) # Phi3^2 - Phi1^2
            )/(Phi1 + Phi3) - 1/6 * h4^2 * Var4 - 1/2 * h4 * Phi2
end

function I_REAL(P, h1, h4)

    P = BigFloat(P)
    h1 = BigFloat(h1)
    h4 = BigFloat(h4)

    hh = h1^2 + h4^2
    R = sqrt(P^2 + hh)
    R1 = sqrt(P^2 + h1^2)
    hh1 = h4^2

    Phi1 = 0.5 * log((P + R)^2 / hh) / P
    Phi2 = atan(h1 * P / (hh + h4 * R)) / (h1 * P)
    Phi3 = 0.5 / R1 * log((R1 + R)^2 / hh1)

    return ((h1^2 - 3 * h4^2) * Phi1 - h4 * (3 * h1^2 - h4^2) * Phi2+ 3 * h4^2 * Phi3 - h4^2 / (R + h4)) / (6 * h1^2)
end
#=
h1 = BigFloat(h1)
h4 = BigFloat(h4)
P = BigFloat(P)
=#
a = (sqrt(P^2+h1^2)+sqrt(P^2+h1^2+h4^2))^2/h4^2 # in Phi3
b = (P+sqrt(P^2+h1^2+h4^2))^2/(h1^2+h4^2) # in Phi1

function new(P, h1, h4, n, m) # n Sqrt, m Ln

    hh = h1^2 + h4^2
    R = sqrt(P^2 + hh)
    R1 = sqrt(P^2 + h1^2)
    RR = P^2 + hh
    hh1 = h4^2

    Phi1 = 0.5 * log((P + R)^2 / hh) / P
    Phi2 = atan(h1 * P / (hh + h4 * R)) / (h1 * P)
    Phi3 = 0.5 / R1 * log((R1 + R)^2 / hh1)

    a = (R1+R)^2/hh1 # in Phi3
    b = (P+R)^2/hh # in Phi1

    Tayl = taylorSqrt((h1/P)^2, n)
    Ref6 = 2/(h4^2 * (P + R)^2) * (h1^2 * (RR + R * R1) + h4^2 * P * R * (Tayl - 1))
    RefLn = taylorLn(Ref6, m)

    Var4 = 1 / (2 * hh * h4 + h4 * P^2 + (2*h4^2 + h1^2) * R)

    return Phi1/6 + 1/2 * (h4/h1)^2 * (
                (P^2 * log(a*b) * RefLn - h1^2 * log(b)^2) / (4 * P^2 * (P^2 + h1^2)) # Phi3^2 - Phi1^2
            )/(Phi1 + Phi3) - 1/6 * h4^2 * Var4 - 1/2 * h4 * Phi2
end

function Implemented(P, h1, h4, n, m)
    
    hh = h1^2 + h4^2
    hh1 = h4^2
    R = sqrt(P^2 + h1^2 + h4^2)
    RR = P^2 + h1^2 +h4^2
    R1 = sqrt(P^2 + h1^2)

    Phi1 = 0.5 * log((P + R)^2 / hh) / P
    Phi2 = atan(h1 * P / (hh + h4 * R)) / (h1 * P)
    Phi3 = 0.5 / R1 * log((R1 + R)^2 / hh1)

    a = (R1+R)^2/hh1 # in Phi3
    b = (P+R)^2/hh # in Phi1

    #taylor series sqrt(1+x)
    sum = 0
    for i in 1:n
        sum += 2 / i * binomial(2 * (i - 1), i - 1) * (-(h1/P)^2 / 4)^i
    end
    TaylSqrt = sum

    Ref6 = 2/(h4^2 * (P + R)^2) * (h1^2 * (RR + R * R1) - h4^2 * P * R * TaylSqrt)
    
    #taylor series ln(1+x)
    sum = 0
    for i in 1:m
        sum -= (-1)^i * (Ref6)^i / i
    end
    RefLn = sum

    Var4 = 1 / (2 * hh * h4 + h4 * P^2 + (2*h4^2 + h1^2) * R)

    return Phi1/6 + 1/2 * (h4/h1)^2 * (
        (P^2 * log(a*b) * RefLn - h1^2 * log(b)^2) / (4 * P^2 * (P^2 + h1^2)) # Phi3^2 - Phi1^2
    )/(Phi1 + Phi3) - 1/6 * h4^2 * Var4 - 1/2 * h4 * Phi2
end

for j in 0:50
    h1 = 2^j * 1e-12
    h4 = 100
    P = 0.01

    RelErr = abs(Implemented(P, h1, h4, 5, 3) - I_REAL(P, h1, h4)) / I_REAL(P, h1, h4)
    println(h1, ", ", RelErr)
end