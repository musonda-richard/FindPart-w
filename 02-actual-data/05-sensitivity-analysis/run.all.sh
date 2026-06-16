#!/bin/bash

set -e

mus="2.0 3.0 4.0 5.0 6.0"
vars="1.0 2.19 3.0" 

for mu in $mus
do
    for var in $vars
    do
        outfile="rgs-AIC_mu${mu}_var${var}.csv"

        if [ -s "$outfile" ]; then
            echo "Skipping mu=${mu}, var=${var} because $outfile already exists"
            continue
        fi

        alpha=$(awk -v mu="$mu" -v var="$var" 'BEGIN {print (mu * mu) / var}')
        theta=$(awk -v mu="$mu" -v var="$var" 'BEGIN {print var / mu}')

        echo "Running mu=${mu}, var=${var}, alpha=${alpha}, theta=${theta}"

        make clean
        make ALPHA=$alpha THETA=$theta

        if [ ! -s rgs-AIC.csv ]; then
            echo "ERROR: rgs-AIC.csv was not created for mu=${mu}, var=${var}"
            exit 1
        fi

        mv rgs-AIC.csv "$outfile"

    done
done

echo "file,mu,var,alpha,theta,rgs,sha1,ll,num_pars,AIC,lb_agg_AIC,has_chance" > combined_results.csv

for f in rgs-AIC_mu*_var*.csv
do

    if [ ! -f "$f" ]; then
        continue
    fi

    base=$(basename "$f" .csv)

    mu=$(echo "$base" | cut -d'_' -f2 | sed 's/mu//')
    var=$(echo "$base" | cut -d'_' -f3 | sed 's/var//')

    if [ -z "$mu" ] || [ -z "$var" ]; then
        echo "Skipping $f because mu or var could not be parsed"
        continue
    fi

    alpha=$(awk -v mu="$mu" -v var="$var" 'BEGIN {print (mu * mu) / var}')
    theta=$(awk -v mu="$mu" -v var="$var" 'BEGIN {print var / mu}')

    awk -F',' \
        -v file="$f" \
        -v mu="$mu" \
        -v var="$var" \
        -v alpha="$alpha" \
        -v theta="$theta" \
        'NR==2 {print file","mu","var","alpha","theta","$0}' \
        "$f" >> combined_results.csv

done

make clean

echo "Done. Created combined_results.csv"
