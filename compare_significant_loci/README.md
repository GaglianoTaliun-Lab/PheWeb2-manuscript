# don't forget to create a venv and install requirements.txt

# run create_loci_files.sh
sbatch create_loci_files.sh

# create file with all loci across phenotypes for each stratification
cat $(ls loci/*.all.male_loci_list.txt ) > total.all.male_loci_list.txt
cat $(ls loci/*.all.female_loci_list.txt) > total.all.female_loci_list.txt
cat $(ls loci/*.all.combined_loci_list.txt ) > total.all.combined_loci_list.txt

# (Or for european vs all)
cat $(ls loci_combined/*.all.combined_loci_list.txt ) > total.all.combined_loci_list.txt
cat $(ls loci_combined/*.european.combined_loci_list.txt) > total.european.combined_loci_list.txt

# run compare_sig_loci.py
sbatch run_compare_sig_loci.sh
