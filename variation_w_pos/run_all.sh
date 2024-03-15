rm output_variance.txt

bash update_backup.sh
bash populate_variance.sh
python plot.py

cp all_data/0.22/* output_data/
python analysis.py 0.22