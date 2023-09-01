#!/bin/sh

# argument entry:
# 1. number of splines in basis
# 2. degree
# 3. amplitude
# 4. iteration number of iterative non linear method (optional)
# 5. folder name (optional)

j=$(($1-1))
for i in $( eval echo {0..$j} )
do
# do echo $i > Iteration.txt;
    source Go_Greens.sh $1 $2 $3 $i $4 $5
    python Save_Strain_Greens.py $i
# taskset --cpu-list $i source Go.sh

done
# echo 50 > Iteration.txt;
python Make_Greens.py $1 $2 $3 $4 $5

# yes | cp -rf GreenFunc/matrix.txt home/rltam/Documents/RunOutputFiles/Greens(25,3)x5
