data_dir=/Users/zhaoyangjia/Library/CloudStorage/OneDrive-UCSanDiego/Bafna_Lab/packaged_data_for_comparison/Bionano_solve/
output_dir=./omkar_output
exe_dir=/Users/zhaoyangjia/Library/CloudStorage/OneDrive-UCSanDiego/Bafna_Lab/OMKar/
centromere_file=/Users/zhaoyangjia/Library/CloudStorage/OneDrive-UCSanDiego/Bafna_Lab/packaged_data_for_comparison/hg38_centro.txt
cyto_file=${exe_dir}/cytoBand.txt
log_dir=${output_dir}/logs/

mkdir -p ${log_dir}

for folder in ${data_dir}/*
do
	folder_name=$(basename "$folder")
	mkdir -p ${output_dir}/${folder_name}
	echo $folder_name
	python ${exe_dir}/main.py \
		-cnv ${folder}/cnv_calls_exp.txt \
		-smap ${folder}/exp_refineFinal1_merged_filter_inversions.smap \
		-rcmap ${folder}/cnv_rcmap_exp.txt \
		-xmap ${folder}/exp_refineFinal1_merged.xmap \
		-n $folder_name \
		-o ${output_dir}/${folder_name} \
		-centro ${centromere_file} \
		-cyto ${cyto_file} > ${log_dir}/${folder_name}.log.txt
done
