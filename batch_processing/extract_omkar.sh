data_dir=/media/zhaoyang-new/workspace/KarSim/packaged_data_for_comparison/OMKar_output_testbuild/
output_dir=/media/zhaoyang-new/workspace/KarSim/KarComparator/new_data_files/OMKar_testbuild/

for folder in ${data_dir}/*
do
	cd ${folder}
	for txt_file in *.1.txt; do
        # Check if the file is not named *.SV.txt
        if [ "$txt_file" != *SV.txt ]; then
            # Copy the file to the destination folder
            cp "$txt_file" "$output_dir"
            echo "Copied $txt_file to $output_dir"
        fi
    done
done
