include("lib_setpartition.jl")
using DataFrames, CSV, SHA, Dates, Statistics
start_time = now()

if length(ARGS) != 6
    @error "Usage: julia hclust.jl count_file baseline estimates_file ftol_rel alpha theta"
    exit(1)
end

count_file = ARGS[1]
baseline   = ARGS[2]
est_file   = ARGS[3]
FTOL_REL = parse(Float64, ARGS[4])
ALPHA    = parse(Float64, ARGS[5])
THETA    = parse(Float64, ARGS[6])

df = CSV.read(count_file, DataFrame)
variants = names(df)[3:end]

deleteat!(variants, findfirst(isequal(baseline), variants))
insert!(variants, 1, baseline)

const n = length(variants)

println("Number of lineages: ", n)
println("Baseline lineage: ", baseline)

df_est = CSV.read(est_file, DataFrame)

rho = df_est.k

println("Loaded rho estimates")

# INITIAL CLUSTERS
clusters = Vector{Vector{Int64}}()
for i in 1:n
    push!(clusters, [i])
end

results = DataFrame(step = Vector{Int64}(), k = Vector{Int64}(), rgs  = Vector{String}(), sha1 = Vector{String}(), 
    ll = Vector{Float64}(), num_pars  = Vector{Int64}(), AIC = Vector{Float64}())

# CONVERT CLUSTERS -> ASSIGNMENTS
function clusters_to_assignments(clusters, n)
    a = zeros(Int, n)
    for cluster_id in 1:length(clusters)
        members = clusters[cluster_id]
        for m in members
            a[m] = cluster_id
        end
    end
    return a
end

# CENTROID
function centroid(cluster, rho)
    return mean(rho[cluster])
end

# CENTROID DISTANCE
function centroid_distance(c1, c2, rho)
    return abs(centroid(c1, rho) - centroid(c2, rho))
end

# FIND CLOSEST CLUSTERS

function find_closest_clusters(clusters, rho)

    best_i = 0
    best_j = 0

    best_dist = Inf

    k = length(clusters)

    for i in 1:(k-1)
        for j in (i+1):k
            d = centroid_distance(clusters[i], clusters[j], rho)
            if d < best_dist
                best_dist = d
                best_i = i
                best_j = j
            end
        end
    end
    return best_i, best_j, best_dist
end

# EVALUATE PARTITION
function evaluate_partition(clusters, step, count_file, baseline,
        FTOL_REL, ALPHA, THETA,est_file, results)

    a = clusters_to_assignments(clusters, n)
    rgs = string_to_rgs(a)
    prefix = bytes2hex(sha1(join(rgs, "-")))
    filename = prefix * "_loglikelihood.csv"

    if !(isfile(filename))
        run(`julia RelRe.jl -b $baseline -a $ALPHA -t $THETA -i $count_file -d 1 -r $(rgs) -o $(prefix) --ftol_rel $FTOL_REL -v $est_file`)
    end

    vec_s = split(readchomp(`tail -n 1 $filename`), ",")
    ll  = parse(Float64, vec_s[1])
    num_pars = parse(Int64, vec_s[2])
    AIC  = parse(Float64, vec_s[3])

    push!(results,(step,length(clusters), join(rgs, "-"), prefix, ll, num_pars, AIC))

    println("Step = ", step," | k = ", length(clusters)," | AIC = ", round(AIC, digits = 4))
end

# EVALUATE INITIAL PARTITION

println("Evaluating initial partition")

evaluate_partition(clusters,0,count_file,baseline,FTOL_REL,ALPHA,THETA,est_file,results)

# CENTROID HIERARCHICAL CLUSTERING
println("Running centroid hierarchical clustering")

for step in 1:(n - 1)

    i, j, d = find_closest_clusters(clusters, rho)

    println("Merge step ", step, " | distance = ", round(d, digits = 6))

    # merge clusters
    new_cluster = vcat(clusters[i],clusters[j])
    sort!(new_cluster)

    deleteat!(clusters, max(i, j))
    deleteat!(clusters, min(i, j))

    # add merged cluster

    push!(clusters, new_cluster)

    println("Clusters remaining = ", length(clusters))

    # evaluate partition
    evaluate_partition(clusters, step, count_file, baseline, FTOL_REL, ALPHA, THETA, est_file, results)
end

sort!(results, :AIC)

CSV.write("hierarchical-rgs-AIC.csv", results)
println("Finished")
println(now())
