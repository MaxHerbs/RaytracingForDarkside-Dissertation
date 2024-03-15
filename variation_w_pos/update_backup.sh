declare -a arr=("0.00" "0.11" "0.33" "0.44" "0.16" "0.27" "0.22")

for i in "${arr[@]}"
do
    echo Pos: $i
    scp godzilla:/storage/physics/phuddg/intensity_w_position/$i/output_data/*  all_data/$i/
done