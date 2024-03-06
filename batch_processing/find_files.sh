# Define the directories
folder1="/media/zhaoyang-new/workspace/joey_snakemake/virtual_karyotype_snakemake/samples/"
folder2="/media/zhaoyang-new/workspace/KarSim/OMKar_outputs/simulation_final_v3/"
folder3="/media/zhaoyang-new/workspace/KarSim/packaged_data/"

# Iterate through folders in folder1
for folder_name in "$folder1"/*; do
    # Extract the base name of the folder
    base_name=$(basename "$folder_name")

    # Check if the same folder exists in folder2
    if [ -d "$folder2/$base_name" ]; then
        # Copy the folder from folder2 to folder3
        cp -r "$folder2/$base_name" "$folder3/"
        echo "Copied $base_name from folder2 to folder3"
    fi
done
