function Basis_SINDY_PI(x, dx, npoly,ntrig, w)
    nvar = size(x, 2)

    basis = ones(size(x,1), 1)

    if npoly >=1
      basis = hcat(basis, x)
    end

    if npoly >= 2
        for i in 1:nvar
            for j in i:nvar
                basis = hcat(basis, (x[:, i] .* x[:, j]) )
            end
        end
    end


    if npoly >= 3
        for i in 1:nvar
            for j in i:nvar
                for k in j:nvar
                basis = hcat(basis, (x[:,i] .* x[:, j] .* x[:, k]) )
            end
            end
        end
    end


    if npoly >= 4
        for i in 1:nvar
            for j in i:nvar
                for k in j:nvar
                    for l in k:nvar
                basis =  hcat(basis, (x[:,i].* x[:, j] .* x[:,k] .* x[:,l]))
            end
            end
        end
    end
end

if npoly >= 5
    for i in 1:nvar
        for j in i:nvar
            for k in j:nvar
                for l in k:nvar
                    for m in l:nvar
            basis =  hcat(basis, (x[:,i].* x[:, j] .* x[:,k] .* x[:,l] .* x[:,m]))
        end
        end
    end
end
end
end

if npoly >= 6
    for i in 1:nvar
        for j in i:nvar
            for k in j:nvar
                for l in k:nvar
                    for m in l:nvar
                        for n in m:nvar
            basis =  hcat(basis, (x[:,i].* x[:, j] .* x[:,k] .* x[:,l] .* x[:,m] .* x[:,n]))
        end
        end
    end
end
end
end
end

if ntrig >=1
    basis = hcat(basis, sin.(x))
    basis = hcat(basis, cos.(x))
end

if ntrig >=2
    basis = hcat(basis, cos.(x) .* cos.(x))
end

if ntrig >=3
    for i in 1:nvar
        for j in 1:nvar
            basis = hcat(basis, (x[:, i] .* x[:, i] .* sin.(x[:, j]) .* cos.(x[:, j]) ))
        end
    end
end

if ntrig >=4
    for i in 1:nvar
        for j in 1:nvar
            basis = hcat(basis, (x[:, i] .* x[:, i] .* sin.(x[:, j])))
        end
    end

    for i in 1:nvar
        for j in 1:nvar
            basis = hcat(basis, (x[:, i] .* x[:, i] .* cos.(x[:, j])))
        end
    end

    for i in 1:nvar
        for j in 1:nvar
            basis = hcat(basis, (sin.(x[:, i]) .* cos.(x[:, j]) ))
        end
    end
end

    pin = size(basis,2)

    basis_noderivative = deepcopy(basis[:, 1:pin])

    basis_deep = zeros(size(basis,1), 2*size(basis,2))

    basis_deep[:, 1:size(basis,2)] = basis_noderivative

    for k=1:pin
        basis_deep[:,pin + k] = dx[:, w].* basis_noderivative[:, k]
    end

    basis_right = zeros(size(basis_deep, 1), size(basis_deep, 2)-1)

    for i = 1:size(basis_noderivative,2)
        basis_right[:, i] = basis_deep[:,i]
    end

    for i = size(basis_noderivative,2)+1:size(basis_deep, 2)-1
        basis_right[:, i] = basis_deep[:,i+1]
    end

    basis_guess = dx[:, w]

    return (basis_right, basis_guess)

end
