# argument entry:
# 1. output folder name
# 2. number of splines in basis
# 3. degree
# 4. amplitude

source Go_ISO.sh $1 $2 $3 $4
python Inversion.py $1 $2 $3 $4

for n in {1..10}
do
    source Greens.sh $2 $3 $4 $n $1
    python Inversion.py $1 $2 $3 $4 $n
done

# one last file to make plots?

source Go_Last.sh $1 $2 $3 $4
python Plots_Iterative.py $1 $2 $3 $4 $n