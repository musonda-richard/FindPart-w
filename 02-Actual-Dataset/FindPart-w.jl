include("lib_setpartition.jl")
using DataFrames, CSV, SHA, Base.Threads, Dates
start_time = now()

if length(ARGS) != 8 
    @error "Usage: julia explore_beam.jl count_file baseline estimates_file loglikelihood_file"
    exit(1)
end
count_file = ARGS[1]
baseline = ARGS[2]
est_file = ARGS[3]
ll_file = ARGS[4]
FTOL_REL = parse(Float64, ARGS[5])
ALPHA = parse(Float64, ARGS[6])
THETA = parse(Float64, ARGS[7])
width = parse(Int64, ARGS[8])

df = CSV.read(count_file, DataFrame)
variants = names(df)[3:end]
deleteat!(variants, findfirst(isequal(baseline), variants))
insert!(variants, 1, baseline)
@show variants

# n is the number of variants, i.e number of subjects + 1 
const n = length(variants)

function separate_rgs(vec_rgs::Vector{Vector{Int64}})
    #run RelRe Program for each rgs in vec_rgs and save rgs and its AIC values to df_AIC
    num = length(vec_rgs) 
    println("separating $num rgs")
    println(now())
    flush(stdout)
    lock = Threads.SpinLock()
    Threads.@threads :greedy for i in 1:num
        rgs = vec_rgs[i]
        tid = Threads.threadid()
        println(string(now()) * " Index i=$i was assigned to thread $tid")
        flush(stdout)
        prefix = bytes2hex(sha1(join(rgs, "-")))
        filename = prefix * "_loglikelihood.csv"
        log = prefix * ".log"
        if !isfile(filename)
            run(pipeline(`julia RelRe.jl -b $baseline -a $ALPHA -t $THETA -i $count_file -d 1 -u -r $rgs -o $prefix --ftol_rel $FTOL_REL -v $est_file`,log))
        end
        println(string(now()) * " Calculation for i=$i was done by thread $tid")
        flush(stdout)
        vec_s = split(readchomp(`tail -n 1 $filename`), ",")
        ll = parse(Float64,vec_s[1])
        num_pars = parse(Int64,vec_s[2])
        AIC = parse(Float64,vec_s[3])
        #rm(prefix * "_loglikelihood.csv")
        #rm(prefix * "_estimates.csv")
        #rm(prefix * "_blocks.csv")
        #the lower bound of AIC of its aggregations
        lb_agg_AIC = AIC - 2 * maximum(rgs) # the most important formula in this program
        Base.Threads.lock(lock) do
            #Add rgs to vec_rgs_no_chance if the rgs has no chance or to vec_rgs_to_check otherwise
            has_chance = true
            if lb_agg_AIC > AIC_min
                push!(vec_rgs_no_chance, rgs)
                has_chance = false
            else
                push!(vec_rgs_to_check, rgs)
                push!(vec_AIC_to_check, AIC)
            end
            if AIC < AIC_min
                global AIC_min = AIC
            end
            push!(df_AIC, (join(rgs, "-"),bytes2hex(sha1(join(vec_rgs[i] , "-"))), ll, num_pars, AIC, lb_agg_AIC, has_chance), promote = true)
        end
    end
    if length(vec_AIC_to_check) != 0
        w = minimum([width, length(vec_AIC_to_check)])
        AIC_threshold = sort(vec_AIC_to_check)[w]
        indices_top_w = findall(x -> x <= AIC_threshold, vec_AIC_to_check) 
        global vec_rgs_to_check = vec_rgs_to_check[indices_top_w]
    end
end

function aggregations_survived()
    vec_survived = Vector{Vector{Int64}}()
    # add all aggregations of rgs in vec_rgs_to_check to vec_rgs
    for i in 1:length(vec_rgs_to_check)
        result = aggregations_of_blocks(vec_rgs_to_check[i])
        lock = Threads.SpinLock()
        Threads.@threads :static for j in 1:length(result)
            tid = Threads.threadid()
            chance = true
            for k in 1:length(vec_rgs_no_chance)
                if is_an_aggregation_of(result[j], vec_rgs_no_chance[k])
                    chance = false
                    break
                end
            end
            if chance
                Base.Threads.lock(lock) do
                    push!(vec_survived, result[j])
                end
            end
        end
    end
    vec_survived = unique(vec_survived)
    return vec_survived
end

vec_rgs_no_chance = Vector{Vector{Int64}}()
vec_rgs_to_check = Vector{Vector{Int64}}()
df_AIC = DataFrame(rgs = Vector{String}(),sha1 = Vector{String}(),ll = Vector{Float64}(), num_pars = Vector{Int64}(), AIC = Vector{Float64}(), lb_agg_AIC = Vector{Float64}(), has_chance = Vector{Bool}())
#add new column named with has_chance::Bool

println("reading level 1")
flush(stdout)
rgs_finest = collect(0:n-1)
df_ll = CSV.read(ll_file, DataFrame)
AIC_finest = df_ll.AIC[1]
ll_finest = df_ll.maxll[1]
num_pars_finest = df_ll.num_pars[1]
has_chance = true
lb_agg_AIC_finest = AIC_finest - 2 * maximum(rgs_finest) #Lower bound of the rgs_finest aggregations
AIC_min = AIC_finest

# store the AIC and lower bound AIC of its aggregations 
push!(df_AIC, (join(rgs_finest, "-"),bytes2hex(sha1(join(rgs_finest, "-"))), ll_finest, num_pars_finest, AIC_finest, lb_agg_AIC_finest, has_chance), promote = true)

# Main loop from second level to final level 
for level in 2:n
    println("\n" * "analysing level " * string(level))
    println(now())
    flush(stdout)
    if level == 2
        println("creating aggregations")
        println(now())
        flush(stdout)
        global vec_rgs = aggregations_of_blocks(rgs_finest)
    else
        println("creating aggregations")
        println(now())
        flush(stdout)
        global vec_rgs = aggregations_survived()
    end
    println(string(length(vec_rgs)) * " RGSs will be analysed")
    println(now())
    flush(stdout)
    if length(vec_rgs) == 0
        println("no rgs to analyse")
        flush(stdout)
        break
    end
    global vec_rgs_to_check = Vector{Vector{Int64}}()
    global vec_AIC_to_check = Vector{Float64}()
    separate_rgs(vec_rgs)
    println("the number of survived RGSs: " * string(length(vec_rgs_to_check)))
    println(now())
    flush(stdout)
end
sort!(df_AIC, :AIC)
CSV.write("rgs-AIC.csv",df_AIC)
df_no_chance = DataFrame(rgs = vec_rgs_no_chance)
df_no_chance.rgs = join.(df_no_chance.rgs, "-")
df_no_chance.sha1 = [bytes2hex(sha1(rg)) for rg in df_no_chance.rgs]
CSV.write("rgs_no_chance.csv", df_no_chance)
println("csv files were created")
println(now())
flush(stdout)

elapsed = Dates.value(now() - start_time)
df_elapsed = DataFrame(milliseconds=elapsed)
CSV.write("elapsed.csv",df_elapsed)
