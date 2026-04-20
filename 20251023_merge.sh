#!/bin/bash

# Create the output directory if it doesn't exist
mkdir -p 04_Merged

# Loop through all ReadsPerGene.out.tab files in 03_Alignment/
for file in 03_Alignment/*_ReadsPerGene.out.tab; do
    filename=$(basename "$file")  # Extract just the filename
    samplename=$(echo "$filename" | sed 's/_ReadsPerGene.out.tab//')

    # Extract gene IDs and read counts (column 4) and save to a new file
    cut -f 1,4 "$file" > 04_Merged/"$samplename"-extractedReadCounts.tab

    # Add header with sample name
    sed -i '1s/^/gene_ID\t'"$samplename"'\n/' 04_Merged/"$samplename"-extractedReadCounts.tab
done

echo "Extraction complete. Output stored in 04_Merged directory."
