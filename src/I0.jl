function I0(P, h1, h2, h3, h4)
    #calculate 1D PBF

    h1 = BigFloat(h1)
    h2 = BigFloat(h2)
    h3 = BigFloat(h3)
    h4 = BigFloat(h4)
    P = BigFloat(P)

    hh = h1^2 + h2^2 + h3^2 + h4^2
    R = sqrt(P^2 + hh)
    hh1 = hh - h1^2
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
                    Phi2 = 1 / (hh + h2 * R) # for case 2 (h1 = 0)
                else
                    Phi2 = atan(h1 * P / (hh + h2 * R)) / (h1 * P) # for case 3
                end
                I = (Phi1 - h2 * Phi2) / 6 # Cases 2 & 3
            else
                R1 = sqrt(P^2 + h1^2)
                if h2^2 + h4^2 < zero2 * hh # if h2 = h4 = 0 -> Case 4
                    if h1 < zerotol # catch h1 = 0
                        I = (Phi1 - 2 / (h3 + R)) / 6 # limit of term below as h1 approaches 0
                    else
                        if h1 * P < zerotol
                            Phi2 = 1 / (hh + h3 * R) # catch P = 0
                        else
                            Phi2 = atan(h1 * P / (hh + h3 * R)) / (h1 * P)
                        end
                        Phi3 = 0.5 * hh1 / R1 * log((R1 + R)^2 / hh1)
                        I = ((h1^2 - h3^2) * Phi1 - 2 * h1^2 * h3 * Phi2 + Phi3) / (6 * h1^2) # Case 4
                    end
                elseif h1^2 + h4^2 < zero2 * hh # if h1 = h4 = 0 -> Case 5
                    h = sqrt(hh)
                    R2 = sqrt(P^2 + h2^2)
                    if h3 < zerotol
                        Phi4 = 0 # catch h3 = 0
                    else
                        Phi4 = h3^2 / (h2 * P^2) * ((R2 / h2 * log((R2 + R) / h3) - log((h2 + h) / h3)))
                    end
                    I = (hh / h2^2 * Phi1 - 1 / (R + h) - Phi4) / 6 # Case 5
                # At this stage h4 =/= 0 and h3 = 0.
                # Pick the last case so that one doesn't divide by 0.
                elseif h1 > 0 # if h2 = h3 = 0 -> Case 6
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
                elseif h2 > 0 # if h1 = h3 = 0 -> Case 7
                    h = sqrt(hh)
                    R2 = sqrt(P^2 + h2^2)
                    Phi2 = atan(h2 * P / (hh + h4 * R)) / P
                    Phi4 = h4^2 / (h2 * P^2) * ((R2 / h2 * log((R2 + R) / h4) - log((h2 + h) / h4)))
                    I = (
                            (1 + 3 * h4^2 / h2^2) * Phi1 - 2 * (h4 / h2)^3 * Phi2 - 3 * Phi4 +
                            (2 * h4^2 - h2^2) / (h2^2 * (R + h))
                        ) / 6 # Case 7
                else
                    I = 0 # no case detected
                end
            end
        end
    end

    return Float64(I)

end
