using LinearAlgebra
function Sparse_matrix_SINDYPI(Theta, dx, lambda, nvars, maxiter)

normLib = zeros(1, size(Theta,2))

Thetan = deepcopy(Theta)

    for norm_k=1:size(Theta,2)
            normLib[norm_k] = norm(Thetan[:,norm_k]);
            Thetan[:,norm_k] = Thetan[:,norm_k]/normLib[norm_k];
    end


   #Epsilon = Thetan\dx

#Below:Only to be used for ill conditioned matrices, like in the glycolysis example. Check on this
Epsilon = zeros(size(Thetan,2),1)
ldiv!(Epsilon, qr(Thetan), dx)


#=
   Epsilon_c = similar(Epsilon)
   Epsilon_c .= Epsilon
=#

    for i in 1:maxiter
        index_remove = abs.(Epsilon) .<= lambda
        Epsilon[index_remove] .= zero(eltype(Epsilon))

        for j in 1:nvars
            index_retain = @. !index_remove[:,j]
            Epsilon[index_retain, j] = Thetan[:, index_retain] \ dx[:,j]
        end
   #=
    if norm(Epsilon_c - Epsilon, 2) < eps()
        break
    else
        Epsilon_c .= Epsilon
    end
=#
    end

        for norm_k=1:length(Epsilon)
            Epsilon[norm_k,:] = Epsilon[norm_k,:]/normLib[norm_k];
        end

    return Epsilon
end

#=
A. Q factor

qrf.R

Ans = inv(qrf.R)*(qrf.Q)'*dx
=#
