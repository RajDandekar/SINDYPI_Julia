function Table_SINDY_PI(s, epsilon, npoly, ntrig, w)

    nvar = length(s)
    table = []
    #push!(table, "")

    push!(table, "1")

   if npoly >=1
    for i in 1:nvar
        push!(table, s[i])
    end
end

    if npoly >= 2
        for i in 1:nvar
            for j in i:nvar
                push!(table, string(s[i], s[j]))
            end
        end
    end


    if npoly >= 3
        for i in 1:nvar
            for j in i:nvar
                for k in j:nvar
                  push!(table, string(s[i], s[j], s[k]))
            end
            end
        end
    end


    if npoly >= 4
        for i in 1:nvar
            for j in i:nvar
                for k in j:nvar
                    for l in k:nvar
                      push!(table, string(s[i], s[j], s[k], s[l]))

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
                  push!(table, string(s[i], s[j], s[k], s[l], s[m]))

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
                  push!(table, string(s[i], s[j], s[k], s[l], s[m], s[n]))

        end
        end
    end
end
end
end
end

if ntrig >=1
   for i in 1:nvar
       push!(table, string("sin", "", s[i]))
   end

   for i in 1:nvar
       push!(table, string("cos", "", s[i]))
   end
end

if ntrig >=2
   for i in 1:nvar
       push!(table, string("cos^{2}", "", s[i]))
   end
end

if ntrig >=3
    for i in 1:nvar
        for j in 1:nvar
            push!(table, string(s[i], "*", s[i], "sin", "", s[j], "cos", "", s[j]))
        end
    end
end

if ntrig >=4
    for i in 1:nvar
        for j in 1:nvar
            push!(table, string(s[i], "*", s[i], "sin", "", s[j]))
        end
    end

    for i in 1:nvar
        for j in 1:nvar
            push!(table, string(s[i], "*", s[i], "cos", "", s[j]))
        end
    end

    for i in 1:nvar
        for j in 1:nvar
            push!(table, string("sin", "", s[i], "cos", "", s[j]))
        end
    end
end


tablecopy = deepcopy(table)

for i in 1:length(table)-1
    push!(table, string("dot", "(", s[w], ")" , "", tablecopy[i+1]))
end


table = hcat(table, epsilon)


    return table
end
