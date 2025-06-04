function GSorthogonalization_expan(e, a)
    # e is a vector; a is a set of vectors that spans a space
    # the returned values:
    # h is the norm of vector component of e that is perpendicular to span(a); 
    # s is the magnitude of vector components in a such that the component of e that is in span(a) is sum(si0*ai)

    T = promote_type(eltype(e), eltype(a[1]))

    ## Orthogonalization using the Gram-Schmidt procedure.
    # We want to find a set of orthogonal vectors u such that span(u) = span(a)
    # Since vectors a don't have to be linearly independent some u may be 0.
    n = length(a)
    u = Vector{Vector{Float64}}(undef, n)

    u[1] = a[1] # start of process
    
    c = zeros(T, n, n) # coefficients of u[k] to get a[i]
    for i in 2:n
        sum_term = zeros(T, 3) #assumed 3D space, so vector u has 3 elements
        for k in 1:(i - 1)
            ukn = norm(u[k])
            if ukn^2 < zerotol # special treatment of linearly dependent vectors a_j
                c[k, i] = T(0.0)
            else
                c[k, i] = dot(u[k], a[i]) / ukn^2
            end
            sum_term += c[k, i] * u[k]
        end
        u[i] = a[i] - sum_term
    end

    #Find all other projection coefficients c
    for j in 1:n
        for i in 1:j
            ujn = norm(u[j])
            if ujn^2 < zerotol
                c[j, i] = T(0.0)
            else
                c[j, i] = dot(u[j], a[i]) / ujn^2
            end
        end
    end

    #Find sj0
    b = zeros(T, n)
    s = zeros(T, n)
    for j in n:-1:1
        b[j] = dot(u[j], e)
        sum_term = T(0.0)
        Ajj = c[j, j] * norm(u[j])^2
        if Ajj < zerotol
            s[j] = T(0.0)
        elseif j == n #sum_term is 0, so ignore it
            s[j] = b[j] / Ajj
        else
            for i in (j + 1):n #calculate sum_term
                Aji = c[j, i] * norm(u[j])^2
                sum_term += Aji * s[i]
            end
            s[j] = (b[j] - sum_term) / Ajj
        end
    end

    #Find h
    sum_term = zeros(T, 3) #same as above, 3D space so vector u has 3 elements
    for i in 1:n
        sum_term += s[i] * a[i]
    end
    h = e - sum_term
    h = norm(h)

    return (h, s)
end
