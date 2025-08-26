function I0(P, h1, h2, h3, h4)
    #calculate 1D PBF

    hh = h1^2 + h2^2 + h3^2 + h4^2
    R = sqrt(P^2 + hh)
    hh1 = h2^2 + h3^2 + h4^2
    zero2 = zerotol * zerotol

    if hh < zero2 * P^2 # if h1 = h2 = h3 = h4 = 0 -> Singularity
        Phi1 = log(P) / P
        I = Phi1 / 6
    else
        Phi1 = 0.5 * log((P + R)^2 / hh) / P
        if hh1 < zero2 * hh # if h2 = h3 = h4 = 0 -> Case 1
            I = Phi1 / 6 # Case 1
        else
            if h3^2 + h4^2 < zero2 * hh # if h3 = h4 = 0 -> Case 2 or 3
                if h1 * P == 0
                    Phi2 = 1/(hh + h2 * R) # for case 2 (h1 = 0)
                else
                    Phi2 = atan(h1 * P / (hh + h2 * R)) / (h1 * P) # for case 3
                end
                I = (Phi1 - h2 * Phi2) / 6 # Cases 2 & 3
            else
                R1 = sqrt(P^2 + h1^2)
                if h2^2 + h4^2 < zero2 * hh # if h2 = h4 = 0 -> Case 4
                    if h1 < 1e-4 * h3 # catch h1 = 0
                        RR = P^2 + h1^2 + h3^2
                        RR1 = P^2 + h1^2

                        I = (
                            1/P * log((P + R)/h3) - 2/(R + h3)
                            + h3^2/2 * ((1/R - 1/R1 * log((R1 + R)/h3)) / RR1 - 1/(RR*P + P^2*R) + 1/(P*hh))
                            ) / 6
                    else
                        if h1 * P < zerotol
                            Phi2 = 1/(hh + h3 * R) # catch P = 0
                        else
                            Phi2 = atan(h1 * P / (hh + h3 * R)) / (h1 * P)
                        end
                        Phi3 = 0.5 * hh1/R1 * log((R1 + R)^2/hh1)

                        I = ((h1^2 - h3^2) * Phi1 - 2 * h1^2 * h3 * Phi2 + Phi3) / (6 * h1^2) # Case 4
                    end
                elseif h1^2 + h4^2 < zero2 * hh # if h1 = h4 = 0 -> Case 5
                if h2 < 1e-4 * h3 # catch h2 = 0
                    R3 = sqrt(P^2 + h3^2)
                    # falsch lol
                    I = (
                        1/P * log((P + R3)/h3) - 1/(R3 + h3)
                        - 1/2 * (h3/P)^2 * (2*h3 + 1/(P + R3) + 1/P * log((P + R3)/h3) + P/h3^2 * (P + R3))
                        ) / 6
                else
                    h = sqrt(hh)
                    R2 = sqrt(P^2 + h2^2)
                    if h3 < zerotol
                        Phi4 = 0 # catch h3 = 0
                    else
                        Phi4 = 1/P^2 * (h3/h2)^2 * ((R2 * log((R2 + R)/h3) - h2 * log((h2 + h)/h3)))
                    end
                    I = ((1 + (h3/h2)^2) * Phi1 - 1/(R + h) - Phi4) / 6 # Case 5
                end
                # At this stage h4 =/= 0 and h3 = 0.
                # Pick the last case so that one doesn't divide by 0.
                elseif h1 > 0 # if h2 = h3 = 0 -> Case 6
                    if h1 * P < zerotol
                        Phi2 = 1/(hh + h4 * R) # catch P = 0
                    else
                        Phi2 = atan(h1 * P / (hh + h4 * R)) / (h1 * P)
                    end
                    Phi3 = 0.5 / R1 * log((R1 + R)^2 / hh1)

                    RR = P^2 + h1^4 + h4^2

                    if h1 < 1e-5 * h4 #&& h1 < P
                        R4 = sqrt(P^2 + h4^2)

                        I = 1/6 * log((P + R4)/h4) / P + 1/4 * (h4/P)^2 * (
                            1/(P + R4) - 1/P * log((P + R4)/h4) + P/h4^2
                        ) - 1/18 * 1/(h4 + R4)^3 * (3 * (h4^2 + h4 * R4) + P) - 1/2 * 1/(h4 + R4)
                    elseif false #h1 < 0.1 * P
                        Phi2 = atan(h1 * P / (hh + h4 * R)) / (h1 * P)

                        a = (R1+R)^2/hh1 # in Phi3
                        b = (P+R)^2/hh # in Phi1

                        #taylor series sqrt(1+x)
                        sum = 0
                        for i in 1:10
                            sum += 2 / i * binomial(2 * (i - 1), i - 1) * (-(h1/P)^2 / 4)^i
                        end
                        TaylSqrt = sum

                        Ref6 = 2/(h4^2 * (P + R)^2) * (h1^2 * (RR + R * R1) - h4^2 * P * R * TaylSqrt)
                        
                        #taylor series ln(1+x)
                        sum = 0
                        for i in 1:5
                            sum -= (-1)^i * (Ref6)^i / i
                        end
                        RefLn = sum

                        Var4 = 1 / (2 * hh * h4 + h4 * P^2 + (2*h4^2 + h1^2) * R)

                        I = Phi1/6 + 1/2 * (h4/h1)^2 * (
                            (P^2 * log(a*b) * RefLn - h1^2 * log(b)^2) / (4 * P^2 * (P^2 + h1^2)) # Phi3^2 - Phi1^2
                        )/(Phi1 + Phi3) - 1/6 * h4^2 * Var4 - 1/2 * h4 * Phi2
                    else
                        I = (
                            (h1^2 - 3 * h4^2) * Phi1 - h4 * (3 * h1^2 - h4^2) * Phi2
                            + 3 * h4^2 * Phi3 - h4^2 / (R + h4)
                            ) / (6 * h1^2) # Case 6 (NICHT umformen)
                    end
                elseif h2 > 0 # if h1 = h3 = 0 -> Case 7
                    h = sqrt(h2^2 + h4^2)
                    R2 = sqrt(P^2 + h2^2)
                    Phi2 = atan(h2 * P / (hh + h4 * R)) / P
                    Phi4 = h4^2 / (h2 * P^2) * ((R2 / h2 * log((R2 + R) / h4) - log((h2 + h) / h4)))

                    if h2 < 1e-4 * h4
                        R4 = sqrt(P^2 + h4^2)

                        I = (
                            1/P * log((P + R4)/h4)
                            - 3/2 * h4^2/P^2 * (1/(P + R4) + P/h4^2 + 1/P * log((P + R4)/h4) - 2/h4)
                            + 1/(R4 + h4)^2 * (2/3 * P^2/(R4 + h4) - R4)
                            ) / 6
                    else
                        I = (
                            (1 + 3 * h4^2 / h2^2) * Phi1 - 2 * (h4 / h2)^3 * Phi2 - 3 * Phi4 +
                            (2 * h4^2 - h2^2) / (h2^2 * (R + h))
                            ) / 6 # Case 7
                    end
                else
                    I = 0 # no case detected
                end
            end
        end
    end

    #debugging
    #println("I = ", I)

    return I

end
