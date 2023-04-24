#!/bin/bash

cd release
cmake -DCMAKE_BUILD_TYPE=Release ..
make

t_values=(0 0.4 0.53 0.62 0.7 0.76 0.82 0.88 0.94 0.9999)
# alpha = 1/t (inf <=> t=0)
alpha_values=(inf 2.5 1.89 1.61 1.43 1.32 1.22 1.14 1.06 1.0001)
ple_values=(2.1 2.25 2.4 2.55 2.75 3 3.35 3.9 5.1 25)


for i in 0 1 2 3 4 5 6 7 8 9
do
  for j in 0 1 2 3 4 5 6 7 8 9
  do
    for k in 1 2 3 4 5
    do
      ./gensatgirg -n 50000 -m 250000 -ple ${ple_values[i]} -t ${t_values[j]} -wseed $((k*10)) -ncseed $((k*110)) -cseed $((k*1000)) -eseed $((k*10000)) -file "ple${i}_t${j}_satgirg_${k}" -dot 0 -edge 1
      ./gengirg -n 50000 -ple ${ple_values[i]} -alpha ${alpha_values[j]} -deg 10 -wseed $((k*10)) -pseed $((k*110)) -sseed $((k*10000)) -file "ple${i}_t${j}_girg_${k}" -dot 0 -edge 1
    done
  done
done