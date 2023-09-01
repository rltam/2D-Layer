# argument entry:
# 1. output folder name
# 2. number of splines in basis
# 3. degree
# 4. amplitude

source Go_ISO.sh $1 $2 $3 $4
# python Inversion_Offset.py $1 $2 $3 $4
python Inversion.py $1 $2 $3 $4