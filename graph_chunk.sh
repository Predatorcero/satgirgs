#!/bin/bash

#SBATCH -A nicola.koessler
#SBATCH --partition sorcery
#SBATCH -C GPU_MEM:80GB
#SBATCH --gpus=1

cd release
cmake -DCMAKE_BUILD_TYPE=Release ..
make

n=50000
m=250000
avg_degree=10
d=2
t_values=(0 0.4 0.53 0.62 0.7 0.76 0.82 0.88 0.94 0.9999)
# alpha = 1/t (inf <=> t=0)
alpha_values=(inf 2.5 1.89 1.61 1.43 1.32 1.22 1.14 1.06 1.0001)
ple_values=(2.1 2.25 2.4 2.55 2.75 3 3.35 3.9 5.1 25)

declare -a combinations
index=0

for i in `seq 0 9`
do
  for j in `seq 0 9`
  do
    for k in `seq 1 5`
    do
      first_seed=$RANDOM
      second_seed=$RANDOM
      third_seed=$RANDOM
      fourth_seed=$RANDOM
      fifth_seed=$RANDOM
      sixth_seed=$RANDOM
      seventh_seed=$RANDOM
      #./gensatgirg -n ${n} -m ${m} -ple ${ple_values[i]} -t ${t_values[j]} -wseed ${first_seed} -ncseed ${second_seed} -cseed ${third_seed} -eseed ${fourth_seed} -file "n${n}_m${m}_ple${i}_t${j}_dimensions${d}_wseed${first_seed}_ncseed${second_seed}_cseed${third_seed}_eseed${fourth_seed}_satgirg${k}" -dot 1 -edge 1
      #./gengirg -n ${n} -ple ${ple_values[i]} -alpha ${alpha_values[j]} -deg ${avg_degree} -wseed ${fifth_seed} -pseed ${sixth_seed} -sseed ${seventh_seed} -file "n${n}_deg${avg_degree}_ple${i}_t${j}_dimensions${d}_wseed${fifth_seed}_pseed${sixth_seed}_sseed${seventh_seed}_girg${k}" -dot 1 -edge 1

      combinations[$index]="$i $j $k $first_seed $second_seed $third_seed $fourth_seed $fifth_seed $sixth_seed $seventh_seed"
      index=$((index + 1))
    done
  done
done

parameters=(${combinations[${SLURM_ARRAY_TASK_ID}]})

i=${parameters[0]}
j=${parameters[1]}
k=${parameters[2]}
first_seed=${parameters[3]}
second_seed=${parameters[4]}
third_seed=${parameters[5]}
fourth_seed=${parameters[6]}
fifth_seed=${parameters[7]}
sixth_seed=${parameters[8]}
seventh_seed=${parameters[9]}

./gensatgirg -n ${n} -m ${m} -ple ${ple_values[i]} -t ${t_values[j]} -wseed ${first_seed} -ncseed ${second_seed} -cseed ${third_seed} -eseed ${fourth_seed} -file "n${n}_m${m}_ple${i}_t${j}_dimensions${d}_wseed${first_seed}_ncseed${second_seed}_cseed${third_seed}_eseed${fourth_seed}_satgirg${k}" -dot 1 -edge 1
./gengirg -n ${n} -ple ${ple_values[i]} -alpha ${alpha_values[j]} -deg ${avg_degree} -wseed ${fifth_seed} -pseed ${sixth_seed} -sseed ${seventh_seed} -file "n${n}_deg${avg_degree}_ple${i}_t${j}_dimensions${d}_wseed${fifth_seed}_pseed${sixth_seed}_sseed${seventh_seed}_girg${k}" -dot 1 -edge 1


