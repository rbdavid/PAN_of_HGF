#!/bin/bash

# usage:
#	bash submit_multi.sh start_num end_num

prod_id=$(printf "%02d" $1)
echo "Submitting Prod $prod_id"
sed -e "s/aaa/$prod_id/g" < prod_n.sh > temp.sh
bsub temp.sh > temp.txt
cat temp.txt
depend_id=$(python3 grab_jobid.py temp.txt)
echo $depend_id

start_index="$(($1+1))"
end_index=$2

for (( i=$start_index;i<=$end_index;i++ ))
do
	prod_id=$(printf "%02d" $i)
	echo "Submitting Prod $prod_id"
	sed -e "s/aaa/$prod_id/g" < prod_n.sh > temp.sh
	bsub -w "$depend_id" temp.sh > temp.txt
	depend_id=$(python3 grab_jobid.py temp.txt)
done

