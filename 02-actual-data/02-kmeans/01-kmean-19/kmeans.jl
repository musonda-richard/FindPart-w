include("lib_setpartition.jl")
using DataFrames, CSV, Clustering, SHA, Base.Threads, Dates
start_time = now()

if length(ARGS) != 8
    @error "Usage: julia kmeans_threads.jl count_file baseline estimates_file loglikelihood_file"
    exit(1)
end

count_file = ARGS[1]
baseline = ARGS[2]
est_file = ARGS[3]
ll_file = ARGS[4]
FTOL_REL = parse(Float64, ARGS[5])
ALPHA = parse(Float64, ARGS[6])
THETA = parse(Float64, ARGS[7])
num_trials = parse(Int64, ARGS[8])

df = CSV.read(count_file, DataFrame)
variants = names(df)[3:end]
deleteat!(variants, findfirst(isequal(baseline), variants))
insert!(variants, 1, baseline)
@show variants

const n = length(variants)

# save lineage order for verification
CSV.write("variant_order.csv", DataFrame(Index = 1:n, Lineage = variants))

df_est = CSV.read(est_file, DataFrame)
k_values = df_est.k
k_matrix = hcat(k_values)
k_matrix_row = reshape(k_matrix, 1, n)

df_AIC = DataFrame(k = Vector{Int64}(), rgs = Vector{String}(), sha1 = Vector{String}(), ll = Vector{Int64}(), num_pars = Vector{Int64}(), AIC = Vector{Float64}())

lock = Threads.SpinLock()

# store assignments
assignments_dict = Dict{Int64, Vector{Vector{Int}}}()

Threads.@threads for k in 1:n
    println("k = " * string(k))
    flush(stdout)

    Base.Threads.lock(lock)
    try
        if k == 1
            for i in 1:num_trials
                rgs = fill(0,n)
                prefix = bytes2hex(sha1(join(rgs, "-")))
                filename = prefix * "_loglikelihood.csv"
                if !(isfile(filename))
                    run(pipeline(`julia RelRe.jl -b $baseline -a $ALPHA -t $THETA -i $count_file -d 1 -r $(rgs) -o $(prefix) --ftol_rel $FTOL_REL -v $(est_file)`, devnull))
                end
                vec_s = split(readchomp(`tail -n 1 $filename`), ",")
                ll = parse(Float64,vec_s[1])
                num_pars = parse(Int64,vec_s[2])
                AIC = parse(Float64, vec_s[3])
                push!(df_AIC, (k, join(rgs, "-"),bytes2hex(sha1(join(rgs, "-"))), ll, num_pars, AIC), promote = true)
            end

        elseif k < n
            for i in 1:num_trials
                r = kmeans(k_matrix_row, k)
                a = assignments(r)

                # save assignments
                if !haskey(assignments_dict, k)
                    assignments_dict[k] = Vector{Vector{Int}}()
                end
                push!(assignments_dict[k], a)

                rgs = string_to_rgs(a)
                prefix = bytes2hex(sha1(join(rgs, "-")))
                filename = prefix * "_loglikelihood.csv"
                if !(isfile(filename))
                    run(pipeline(`julia RelRe.jl -b $baseline -a $ALPHA -t $THETA -i $count_file -d 1 -r $(rgs) -o $(prefix) --ftol_rel $FTOL_REL -v $est_file`, devnull))
                end
                vec_s = split(readchomp(`tail -n 1 $filename`), ",")
                ll = parse(Float64,vec_s[1])
                num_pars = parse(Int64,vec_s[2])
                AIC = parse(Float64, vec_s[3])
                push!(df_AIC, (k, join(rgs, "-"),bytes2hex(sha1(join(rgs, "-"))), ll, num_pars, AIC), promote = true)
            end

        else # k == n
            for i in 1:num_trials
                rgs = collect(0:n-1)
                prefix = bytes2hex(sha1(join(rgs, "-")))
                filename = prefix * "_loglikelihood.csv"
                if !(isfile(filename))
                    run(pipeline(`julia RelRe.jl -b $baseline -a $ALPHA -t $THETA -i $count_file -d 1 -r $(rgs) -o $(prefix) --ftol_rel $FTOL_REL -v $est_file`, devnull))
                end
                vec_s = split(readchomp(`tail -n 1 $filename`), ",")
                ll = parse(Float64,vec_s[1])
                num_pars = parse(Int64,vec_s[2])
                AIC = parse(Float64, vec_s[3])
                push!(df_AIC, (k, join(rgs, "-"),bytes2hex(sha1(join(rgs, "-"))), ll, num_pars, AIC), promote = true)
            end
        end
    finally
        Base.Threads.unlock(lock)
    end
end

sort!(df_AIC, :AIC)
CSV.write("kmeans-rgs-AIC.csv", df_AIC)

# consensus function
function consensus_matrix(assignments_list, n)
    M = zeros(Float64, n, n)
    n_runs = length(assignments_list)

    for a in assignments_list
        for i in 1:n
            for j in 1:n
                if a[i] == a[j]
                    M[i, j] += 1
                end
            end
        end
    end

    return M ./ n_runs
end

# WRITE CONSENSUS WITH LINEAGE NAMES
println("Available k values: ", collect(keys(assignments_dict)))

for k_selected in sort(collect(keys(assignments_dict)))  
    println("Computing consensus for k = ", k_selected)

    C = consensus_matrix(assignments_dict[k_selected], n)

    dfC = DataFrame(C, Symbol.(variants))

    # add row labels
    dfC = hcat(DataFrame(Lineage = variants), dfC)

    CSV.write("consensus_k$(k_selected).csv", dfC)
end

println(now())
flush(stdout)
