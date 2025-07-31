using CSV, DataFrames
# Get values from bash
filename = ARGS[1]
df_params = DataFrame(CSV.File(filename))
values = df_params.k

# Perform unique value mapping
unique_values = unique(values)
value_to_int = Dict(value => i for (i, value) in enumerate(unique_values))
mapped_values = [value_to_int[v] for v in values]

function string_to_rgs(vec::Vector{Int64})
    result = copy(vec)

    if result[1] != 0
        new_max = maximum(result) + 1
        index = findall(x -> x == 0, result)
        result[index] .= new_max
        index = findall(x -> x == result[1], result)
        result[index] .= 0
    end

    for i in 2:length(result)
        if result[i] <= (maximum(result[1:(i-1)]) + 1) && result[i] >= 0
            continue
        else
            num_to_be_used = maximum(result[1:i-1]) + 1
            index = findall(x -> x == num_to_be_used, result)
            result[index] .= maximum(result) + 1
            index = findall(x -> x == result[i], result)
            result[index] .= num_to_be_used
        end
    end
    return result
end

converted_values = string_to_rgs(mapped_values) 
converted_values_str = join(converted_values, "-") 
println(converted_values_str)

