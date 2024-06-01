build_version=$1

#################MODIFY########################
data_dir=/media/zhaoyang-new/workspace/KarSim/packaged_data_for_comparison/Bionano_solve/
exe_dir=/media/zhaoyang-new/workspace/OMKar/
###############################################

output_dir_prefix=$(pwd)
output_dir=${output_dir_prefix}/builds/${build_version}/
centromere_file=${exe_dir}/hg38_centro.txt
cyto_file=${exe_dir}/cytoBand.txt
log_dir=${output_dir}/logs/
version_register=${output_dir_prefix}/omkar_versions.txt

export PYTHONPATH=../:$PYTHONPATH
mkdir -p ${output_dir}/omkar_output/
mkdir -p $log_dir

cd $exe_dir
commit_code=$(git log -1 --format=%h)
cd "$output_dir_prefix"
echo "${build_version}    ${commit_code}" >> ${version_register}

### Run OMKar
#echo "---------------RUNNING OMKAR-------------------"
#for folder in ${data_dir}/*
#do
#	folder_name=$(basename "$folder")
#	mkdir -p ${output_dir}/omkar_output/${folder_name}/
#	echo $folder_name
#	python ${exe_dir}/main.py \
#		-cnv ${folder}/cnv_calls_exp.txt \
#		-smap ${folder}/exp_refineFinal1_merged_filter_inversions.smap \
#		-rcmap ${folder}/cnv_rcmap_exp.txt \
#		-xmap ${folder}/exp_refineFinal1_merged.xmap \
#		-n ${folder_name} \
#		-o ${output_dir}/omkar_output/${folder_name}/ \
#		-centro ${centromere_file} \
#		-cyto ${cyto_file} > ${log_dir}/${folder_name}.log.txt
#done
#
#### Extract OMKar
#echo "---------------EXTRACTING OMKAR-------------------"
#data_dir=${output_dir}/omkar_output/
#
#mkdir -p ${output_dir}/omkar_paths/
#
#for folder in ${data_dir}/*
#do
#	for txt_file in ${folder}/*.txt; do
#        # Check if the file is not named *.SV.txt
#        if [[ ! "$txt_file" =~ SV\.txt$ && ! "$txt_file" =~ preILP_nodes\.txt$ && ! "$txt_file" =~ preILP_edges\.txt$ ]]; then
#            # Copy the file to the destination folder
#            cp "$txt_file" ${output_dir}/omkar_paths/
#        fi
#    done
#done
#
#### Clustering
#echo "---------------CLUSTERING-------------------"
#python batch_clustering.py \
#  --omkar_output ${output_dir}/omkar_output/ \
#  --omkar_paths_output ${output_dir}/omkar_paths/ \
#  --output_dir ${output_dir}/cluster_files/

### Comparison
echo "---------------COMPARING-------------------"
mkdir -p ${output_dir}/analyses_summary/
python generate_analyses.py \
  --data_dir ${output_dir}/cluster_files/ \
  --output_dir ${output_dir}/analyses_summary/ \
  --omkar_log_dir ${output_dir}/omkar_output/



