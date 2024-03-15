declare -a arr=("0.00" "0.11" "0.33" "0.44" "0.16" "0.27" "0.22")



for i in "${arr[@]}"
do
    cp all_data/$i/* output_data/
    python analysis.py $i
done