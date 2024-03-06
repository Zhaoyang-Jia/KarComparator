executable_dir=/media/zhaoyang-new/workspace/KarSim/KarSimulator/
data_dir=/media/zhaoyang-new/workspace/KarSim/packaged_data/KarSimulator_output_filtered/
output_dir=/media/zhaoyang-new/workspace/KarSim/packaged_data/KarSimulator_output_diploid/

cd $executable_dir

for file in ${data_dir}/*
do
	python Main/KT_haploid_to_diploid.py $file $output_dir
done
