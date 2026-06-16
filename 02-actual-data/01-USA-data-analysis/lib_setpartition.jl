using Combinatorics

# function that checks if supplied vectors follow retricted growth string pattern
function is_restricted_growth(vec::Vector{Int64})
  if vec[1] != 0 
    return false
  end
  for i in 2:length(vec)
    if (vec[i] > 1 + maximum(vec[1:(i - 1)]) || vec[i] < 0)
      return false
    end
  end
  return true
end

function string_to_rgs(vec::Vector{Int64})
  result = Vector{Int64}(vec)
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

# function to convert restricted growth strings into partitions
function blocks_to_partition(vec::Vector{Int64})
  if !is_restricted_growth(vec)
    @error "The vector is not a restricted growth string"
  end
  vec_uniq = unique(vec)
  partition = Vector{Vector{Int64}}(undef, length(vec_uniq))
  for i in 1:length(vec_uniq)
    partition[i] = findall(x -> x == vec_uniq[i], vec)
  end
  return partition
end

# function to check if a given vector is a refinement of another given vector
function is_a_refinement_of(vec_x::Vector{Int64}, vec_y::Vector{Int64})
  if length(vec_x) != length(vec_y)
    @error "Length of the vectors should be the same"
  end
  new_x = blocks_to_partition(vec_x)
  new_y = blocks_to_partition(vec_y)
  for i in 1:length(new_x)
    elem_x = new_x[i]
    positions_y=findall(elem_y->issubset(elem_x,elem_y), new_y)
    if length(positions_y) == 0
      return false
    end
  end
  return true
end

# function to check if a given vector is an aggregation of another given vector
function is_an_aggregation_of(vec_y::Vector{Int64}, vec_x::Vector{Int64})
  return is_a_refinement_of(vec_x, vec_y)
end

# function that converts a given set partition into its respective restricted growth string pattern
function partition_to_blocks(partition::Vector{Vector{Int64}})
  sorted_part = sort_partition(partition)
  max_value = maximum(vcat(sorted_part...))
    
  growth_string = Vector{Int64}(undef, max_value)
  for group in 1:length(sorted_part)
    for element in sorted_part[group]
      growth_string[element] = group - 1
    end
  end

  return growth_string
end

# function that checks if a given vector is a partition of another given vector
function is_partition(X::Vector{Vector{Int64}},V::Vector{Int64})
  # X is a non-empty subset of V such that each element in V is contained in just one element of V
  for x in X
    if length(x) == 0
      return false
    end
  end
  for v in V
    if length(findall(x->v in x,X)) != 1
      return false
    end
  end
  return true
end 

# function that sorts vector of set partitions into ascending order
function sort_partition(X::Vector{Vector{Int64}})
  result = Vector{Vector{Int64}}(undef,length(X))
  for i in 1:length(X)
    result[i] = sort(X[i])
  end
  return sort(result,lt=(a,b) -> a[1]<b[1])
end

# function that unionises all vectors in a vector of partitions into their original universal set
function partition_to_set(X::Vector{Vector{Int64}})
  vec_v = sort(vcat(X...))
  return vec_v
end

# function that enumerates all the possible partitions using restricted growth string partern
# Algorithm V in Stamatelatos 2021
function next_rgs(vec_a::Vector{Int64},vec_b::Vector{Int64})
  n = length(vec_a)
  d = n
  while(vec_a[d] == n-1 || vec_a[d] > vec_b[d])
    d = d - 1
  end
  if(d == 1)
    return false # no more rgs
  else
    vec_a[d] = vec_a[d] + 1
    if(d < n)
      for i in d+1:1:n
        vec_a[i] = 0
        vec_b[i] = maximum([vec_a[i-1], vec_b[i-1]])
      end
    end
  end
  return true # there is another rgs
end

function enumerate_rgs(n::Int64)
  vec_a = fill(0,n)
  vec_b = fill(0,n)
  vec_all = vcat(vec_a)
  while next_rgs(vec_a,vec_b) 
    append!(vec_all,vec_a)
  end
  mat_a = reshape(vec_all, n,:)'
  return mat_a
end

# function that enumerates all the  partitions with at most k groups using restricted growth string partern
# Algorithm W in Stamatelatos 2021
function next_rgs_at_most(vec_a::Vector{Int64},vec_b::Vector{Int64},k::Int64)
  n = length(vec_a)
  d = n
  while(d > 1 && (vec_a[d] == k-1 || vec_a[d] > vec_b[d]))
    d = d - 1
  end
  if(d == 1)
    return false # no more rgs
  else
    vec_a[d] = vec_a[d] + 1
    if(d < n)
      for i in d+1:1:n
        vec_a[i] = 0
        vec_b[i] = maximum([vec_a[i-1], vec_b[i-1]])
      end
    end
  end
  return true # there is another rgs
end

function enumerate_rgs_at_most(n::Int64, k::Int64)
  vec_a = fill(0,n)
  vec_b = fill(0,n)
  vec_all = vcat(vec_a)
  while next_rgs_at_most(vec_a,vec_b,k) 
    append!(vec_all,vec_a)
  end
  mat_a = reshape(vec_all, n,:)'
  return mat_a
end

# function that enumerates only partitions containing exactly k distinct rgs blocks patterns
function next_rgs_exactly_x(vec_a::Vector{Int64},vec_b::Vector{Int64},k::Int64)
  n = length(vec_a)
  result = false
  m = 0
  while m != k-1
    result = next_rgs_at_most(vec_a,vec_b,k)
    m = maximum([vec_a[n],vec_b[n]])
  end
  return result
end

function enumerate_rgs_exactly_x(n::Int64, k::Int64)
  vec_a =vcat(fill(0,n-k), collect(0:(k-1)))
  vec_b = fill(0,n)
  vec_all = vcat(vec_a)
  while next_rgs_exactly_x(vec_a,vec_b,k) 
    append!(vec_all,vec_a)
  end
  mat_a = reshape(vec_all, n,:)'
  return mat_a
end

# faster function that enumerates only partitions containing exactly k distinct rgs blocks patterns
function next_rgs_exactly_y(vec_a::Vector{Int64},vec_b::Vector{Int64},k::Int64)
  n = length(vec_a)
  result = next_rgs_at_most(vec_a,vec_b,k)
  if (maximum([vec_a[n],vec_b[n]]) != k -1)
    for i in n:-1:1
      k0 = k + (i - n - 1)
      if k0 > vec_b[i] && k0 > 0
        vec_a[i] = k0
        vec_b[i] = k0 - 1
      else
        break
      end
    end
  end
  return result
end


function enumerate_rgs_exactly_y(n::Int64, k::Int64)
  vec_a =vcat(fill(0,n-k), collect(0:(k-1)))
  vec_b = fill(0,n)
  vec_all = vcat(vec_a)
  while next_rgs_exactly_y(vec_a,vec_b,k) 
    append!(vec_all,vec_a)
  end
  mat_a = reshape(vec_all, n,:)'
  return mat_a
end

function next_rgs_reverse(vec_a::Vector{Int64})
  n = length(vec_a)
  d = n
  while vec_a[d] == 0 
    d = d - 1
  end
  if d == 1
    return false  # no more rgs
  else
    vec_a[d] = vec_a[d] - 1
    if d < n
      for i in (d+1):n
        vec_a[i] = maximum(vec_a[1:i-1]) + 1
      end
    end
  end
  return true  # there is another rgs
end

function enumerate_rgs_reverse(n::Int64)
  vec_a = collect(0:n-1)
  vec_all = vcat(vec_a)
  while next_rgs_reverse(vec_a)
    append!(vec_all, vec_a)
    if vec_a == fill(0, n)
      break
    end
  end
  mat_a = reshape(vec_all, n,:)'
  return mat_a
end

function enumerate_rgs_reverse_vec(n::Int64)
    vec_rgs = Vector{Vector{Int64}}() 
    vec_a = collect(0:n-1)
    push!(vec_rgs, Vector{Int64}(vec_a))
    while next_rgs_reverse(vec_a)
        push!(vec_rgs, Vector{Int64}(vec_a))
        if vec_a == fill(0, n)
            break
        end
    end
    return vec_rgs
end


# function that returns all possible one step aggregations of a given partition
function aggregations_by_merging(part::Vector{Vector{Int64}})
  result = Vector{Vector{Vector{Int64}}}()
  n = length(part)
  pairs = collect(combinations(1:n,2))
  for pair in pairs
    aggregation = Vector{Vector{Int64}}()
    s = pair[1]
    t = pair[2]
    for i in 1:n
      if i == s
        push!(aggregation, vcat(part[s], part[t]))
      elseif i == t
        continue
      else
        push!(aggregation, part[i])
      end
    end
    push!(result, aggregation)
  end
  return result
end

# function that returns all possible supsets of a given partition
function all_aggregations_of(part::Vector{Vector{Int64}})
  results = Vector{Vector{Vector{Int64}}}()
  push!(results, sort_partition(part))
  for agg_1 in aggregations_by_merging(part)
    for agg_2 in all_aggregations_of(agg_1)
      if agg_2 in results
        continue
      else
        push!(results, agg_2)
      end
    end
  end
  return results
end

# function that returns common aggregations at level n+1 for given two RGSs at level n 
function immediate_common_aggregations(rgs1::Vector{Int64}, rgs2::Vector{Int64})
    @assert is_restricted_growth(rgs1)
    @assert is_restricted_growth(rgs2)
    @assert maximum(rgs1) == maximum(rgs2)
    agg1 = aggregations_of_blocks(rgs1) 
    agg2 = aggregations_of_blocks(rgs2) 
    results = intersect(agg1,agg2) 
    return results 
end

# function that returns the finest common aggregation for given two maximal RGSs
function finest_common_aggregation(rgs1::Vector{Int64}, rgs2::Vector{Int64})
    @assert is_restricted_growth(rgs1)
    @assert is_restricted_growth(rgs2)
    @assert length(rgs1) == length(rgs2)
    s = copy(rgs1)
    t = copy(rgs2)
    l = length(rgs1)
    #@show s t
    #println("")
    while s != t 
        for i in 1:l
            if s[i] == t[i]
                continue
            elseif s[i] < t[i]
                #modify t
                old_t_i = t[i]
                pos_t = findall(x -> x == old_t_i,t)
                t[pos_t] .= s[i]
                pos_t = findall(x -> x > old_t_i,t)
                t[pos_t] = t[pos_t] .- 1
            else # s[i] > t[i]
                #modify s
                old_s_i = s[i]
                pos_s = findall(x -> x == old_s_i,s)
                s[pos_s] .= t[i]
                pos_s = findall(x -> x > old_s_i,s)
                s[pos_s] = s[pos_s] .- 1
            end
            #@show s t
            #println("")
        end
    end
    return s
end

function aggregations_of_blocks(vec::Vector{Int64})
  @assert is_restricted_growth(vec)
  result = Vector{Vector{Int64}}()
  m = maximum(vec)
  pairs = collect(combinations(0:m,2))

  for pair in pairs
    aggregation = Vector{Int64}(vec)
    s = pair[1] # 0
    t = pair[2] # 1
    pos_s = findall(x -> x == s, vec)
    pos_t = findall(x -> x == t, vec)

    for i in 0:m #  2
      if i == s
        continue
      elseif i == t
        index = findall(x -> x == t, aggregation)
        aggregation[index] .= s
      else
        if i > t
          index = findall(x -> x == i, aggregation)
          aggregation[index] .= i - 1
        end
      end
    end
    push!(result, aggregation)
  end
  return result
end
