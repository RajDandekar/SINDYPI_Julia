function SINDYPI(X, DX, Polyorder, Trigorder)
    s = []
    for i = 1:1:size(X,2)
        push!(s, string("x", string(i)))
    end

    number_states_guess = size(X, 2)

    Training_frac = 0.5

    λ_span = (0.1,0.1)
    datasize = 1
    λ = range(λ_span[1], λ_span[2], length = datasize)

    Store_Result = Array{Any}(undef,length(λ), number_states_guess)
    Index_nonzero = Array{Any}(undef,length(λ), number_states_guess)

    Error = zeros(length(λ) , number_states_guess)

    Xi = Array{Any}(undef,length(λ), number_states_guess)

    X_train = X[1:Int(floor(Training_frac*size(X,1))),:]
    DX_train = DX[1:Int(floor(Training_frac*size(DX,1))),:]

    X_test = X[Int(floor(Training_frac*size(X,1)))+1:end,:]
    DX_test = DX[Int(floor(Training_frac*size(DX,1)))+1:end,:]


    for w = 1:number_states_guess
    for i = 1:datasize
    Basis_train = Basis_SINDY_PI(X_train, DX_train, Polyorder[w], Trigorder[w], w)

    Left_train = Basis_train[2]
    Right_train = Basis_train[1]


    Xi[i,w] = Sparse_matrix_SINDYPI(Right_train,Left_train, λ[i], 1, 10)


    String_table= Table_SINDY_PI(s, Xi[i,w], Polyorder[w], Trigorder[w], w)

    Store_Result[i, w] = String_table

    Index_nonzero[i,w] = findall(x->x !=0, Store_Result[i, w][:,2])

    Basis_test =  Basis_SINDY_PI(X_test, DX_test, Polyorder[w], Trigorder[w], w)

    Left_test = Basis_test[2]
    Right_test = Basis_test[1]

    Error[i,w] = norm(Left_test - Right_test[:, Index_nonzero[i,w]] * Xi[i,w][Index_nonzero[i,w]])

    end

    Ind_min = argmin(Error[:, w])

    ODE_Index = Index_nonzero[Ind_min, w]

    ODE_Strings = Store_Result[Ind_min, w][:,1][ODE_Index]

    ODE_values = Xi[Ind_min,w][ODE_Index]

    ODE_Print = []
    push!(ODE_Print, "dot", s[w])
    push!(ODE_Print, "", "=")


    global count = 0
    for i = 1:length(ODE_Strings)
        global count
     if count >0
        if string(round(ODE_values[i]))[1] == '-'
            pop!(ODE_Print)
        end
      end

        push!(ODE_Print, string(string(round(ODE_values[i], digits = 2)), ODE_Strings[i]))
        push!(ODE_Print, "", "+")

        count += 1
    end

    pop!(ODE_Print)

    ODE_Print = join(ODE_Print)

    println("The differential equation for ", s[w], " is"  )

    println(ODE_Print)
    end

end
