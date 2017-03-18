#!/bin/sh

Exe="../../sage"
file="../HFE.sage"
array=( 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 )

for num in "${array[@]}"
do
	average=0
	echo "plaintext length=$num"
	for ((i=10; i >0; i--)) ;
	do
		result=$("$Exe" "$file" "$num")
		average=$(echo "$average + $result"|bc)
	done
	average=$(echo "scale=3; $average/10"|bc) 
	echo "Average time=$average"
done


