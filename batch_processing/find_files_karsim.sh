# Define the directories
folder1="/media/zhaoyang-new/workspace/KarSim/packaged_data/OMKar_output/"
folder2="/media/zhaoyang-new/workspace/KarSim/packaged_data/KarSimulator_output_raw/"
folder3="/media/zhaoyang-new/workspace/KarSim/packaged_data/KarSimulator_output/"

# Iterate through folders in folder1
for folder_name in "$folder1"/*; do
    # Extract the base name of the folder
    base_name=$(basename "$folder_name")
    base_name="${base_name%%.*}"

    echo $base_name
    # Check if the same folder exists in folder2
    if [ -e "$folder2/${base_name}.kt.txt" ]; then
        # Copy the folder from folder2 to folder3
        cp "$folder2/${base_name}.kt.txt" "$folder3/"
        echo "Copied $base_name from folder2 to folder3"
    fi
done
