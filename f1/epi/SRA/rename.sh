#!/bin/bash
#SBATCH -J Rename
#SBATCH -o Rename.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 5  
#SBATCH --mem-per-cpu 2G
#SBATCH --time 01-00:00


# File containing the list of old and new filenames
rename_list=/RAID1/working/R425/lavakau/chromatin/meta.txt

# Check if the rename list file exists
if [[ ! -f $rename_list ]]
then
    echo "Rename list file '$rename_list' not found!"
    exit 1
fi

# Read the rename list and rename files using only the first two columns, skipping the header
cd /RAID1/working/R425/lavakau/chromatin/data1
awk 'NR > 1' "$rename_list" | while read -r line
do
    oldname=$(echo "$line" | awk '{print $1}')
    newname=$(echo "$line" | awk '{print $2}')

    if [[ -f "$oldname" ]]
    then
        if [[ -f "$newname" ]]
        then
            echo "File '$newname' already exists! Skipping rename of '$oldname'."
        else
            mv "$oldname" "$newname"
            echo "Renamed '$oldname' to '$newname'"
        fi
    else
        echo "File '$oldname' not found!"
    fi
done

cd /RAID1/working/R425/lavakau/chromatin/data2
awk 'NR > 1' "$rename_list" | while read -r line
do
    oldname=$(echo "$line" | awk '{print $1}')
    newname=$(echo "$line" | awk '{print $2}')

    if [[ -f "$oldname" ]]
    then
        if [[ -f "$newname" ]]
        then
            echo "File '$newname' already exists! Skipping rename of '$oldname'."
        else
            mv "$oldname" "$newname"
            echo "Renamed '$oldname' to '$newname'"
        fi
    else
        echo "File '$oldname' not found!"
    fi
done