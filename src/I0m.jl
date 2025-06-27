function I0m(P, h1, h2, h3, h4)
    # calculate 1D PBF

    # The double layer integrals can be computed using only the cases where h4 =/= 0.
    # Therefore, only cases 6 and 7 need implementation below.
    # Since h4 is different from 0, it follows that h3 = 0 always. 

    zero2 = zerotol^2
    hh = h1^2 + h2^2 + h4^2
    R = sqrt(P^2 + hh)
    Phi1 = 0.5 * log((P + R) * (P + R) / hh) / P
    R1 = sqrt(P * P + h1 * h1)

    if hh - h4^2 < zero2 * hh # if h1 = h2 = 0
        I = 0.0 # based on the limit of case 6 as h1 approaches 0
    elseif h2 * h2 < zero2 * hh # if h2 = 0 -> Case 6
        if P < zerotol
            Phi2 = h1 / (hh + h4 * R)
        else
            Phi2 = atan(h1 * P / (hh + h4 * R)) / P
        end
        Phi3 = 1 / R1 * log((R1 + R) / h4)
        I = (
            Phi1 + (h1 * h1 - h4 * h4) / (2 * h1 * h4) * Phi2
            - Phi3 + 0.5 / (R + h4)
            ) / (h1 * h1) # Case 6
    elseif h1 * h1 < zero2 * hh # if h1 = 0 -> Case 7
        h = sqrt(hh)
        R2 = sqrt(P * P + h2 * h2)
        Phi2 = atan(h2 * P / (hh + h4 * R)) / P
        Phi4 = h4 * h4 / (h2 * P * P) * ((R2 / h2 * log((R2 + R) / h4) - log((h2 + h) / h4)))
        I = (
            - Phi1 + h4 / h2 * Phi2 + h2 * h2 / (h4 * h4) * Phi4
            - (R2 * R2 / (R + h4) - h2 * h2 / (h + h4)) / (P * P)
            ) / (h2 * h2) # Case 7
    end

    return I
end
