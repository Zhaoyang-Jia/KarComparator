build_version=$1

data_dir=/media/zhaoyang-new/workspace/KarSim/packaged_data_for_comparison/Bionano_solve/
output_dir=/Users/zhaoyangjia/PyCharm_Repos/KarComparator/omkar_analyses_pipeline/builds/${build_version}/
exe_dir=/media/zhaoyang-new/workspace/OMKar/
centromere_file=${exe_dir}/hg38_centro.txt
cyto_file=${exe_dir}/cytoBand.txt
log_dir=${output_dir}/logs/

mkdir -p ${output_dir}/omkar_output/
mkdir -p $log_dir

## Run OMKar
for folder in ${data_dir}/*
do
	folder_name=$(basename "$folder")
	mkdir -p ${output_dir}/omkar_output/${folder_name}
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

## Extract OMKar: update using Linux version
data_dir=${output_dir}/omkar_output/

mkdir -p ${output_dir}/omkar_paths/

for folder in ${data_dir}/*
do
	for txt_file in ${folder}/*.1.txt; do
        # Check if the file is not named *.SV.txt
        if [ "$txt_file" != *SV.txt ]; then
            # Copy the file to the destination folder
            cp "$txt_file" "$output_dir"
            echo "Copied $txt_file to $output_dir"
        fi
    done
done

## Batch_Clustering



