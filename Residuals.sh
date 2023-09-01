# argument entry:
# 1. number of splines in basis
# 2. degree
# 3. amplitude

for x in {0..99};
do
	echo $x
	python MonteCarloCoeffs.py $1
	source Go_ISO.sh $x $1 $2 $3
	python Inversion.py $x $1 $2 $3

done
