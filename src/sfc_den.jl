function build_den_sfc(G::Vector{Vector{Int64}}, num_p::Int64)
    partition_sfc = Dict(i => Vector{Int64}() for i in 1:num_p)
    lookups = Dict(i => -1 for i in 1:length(G))

    num_larger = length(G) % num_p
    size_of_part = length(G) ÷ num_p

    curr_part = 1
    count = 0
    for node in 1:length(G)
        neighs = G[node]
        push!(partition_sfc[curr_part], node)
        lookups[node] = curr_part
        count += 1

        if num_larger > 0
            if (count == (size_of_part + 1) && (curr_part != num_p))
                curr_part += 1
                count = 0
                num_larger -= 1
            end
        else
            if ((count == size_of_part) && (curr_part != num_p))
                curr_part += 1
                count = 0
            end

        end

    end
    return partition_sfc, lookups
end
