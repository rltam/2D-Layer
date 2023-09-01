j=$(($1-1))
for i in $( eval echo {0..$j} )
do
    folder='hello'
    echo $i
    python test.py $j $2 $3

done