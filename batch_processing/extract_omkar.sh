data_dir=omkar_output_temp/
output_dir=../new_data_files/OMKar_testbuild2/

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