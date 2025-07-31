using Distributions,ArgParse ,CSV , DataFrames, Printf
const n = parse(Int64, ARGS[1])
const h = parse(Int64, ARGS[2]) 
const min_k = 0.9
const max_k = 2
const min_qt = 0.01
const max_qt = 0.1
const min_qt_baseline = 0.70
const max_qt_baseline = 0.95
start_date = "2024-01-01"

vec_date = fill(start_date,n)
vec_c = fill(1.0,n)
qt_baseline = rand(Uniform(min_qt_baseline,max_qt_baseline),1)
qt_subjects = rand(Uniform(min_qt,max_qt), n-1)
vec_qt = vcat(qt_baseline, qt_subjects)
vec_qt = vec_qt/sum(vec_qt)

vec_k_temp = vcat(1.0,rand(Uniform(min_k, max_k), h-1))
vec_k = fill(1.0,n) 
while length(unique(vec_k)) != h
    global vec_k = vcat(1.0, sample(vec_k_temp, n-1, replace = true))
end

if length(vec_k) != n
    println("length of vec_k must be $n")
    @assert false
end

variants = map(x -> @sprintf("v%03i", x), 0:(n-1))
df = DataFrame(variant=variants,date=vec_date,c=vec_c,k=vec_k,qt=vec_qt)
CSV.write("parameters.csv",df)
