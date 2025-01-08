#!/bin/bash

# Script to upload ENA reads

# Output file for successful submissions
success_log="add_path_here"

# Clear the previous success log
> "$success_log"

# Iterate through directories
for directory in /path/to/manifest_files/*; do 
    # Iterate through manifest files in each directory
    for manifest in "$directory"/*_manifest.txt; do
        # Capture the full output of the Java command
        output=$(java -jar add_path_to.jar \
            -context reads \
            -manifest "${manifest}" \
            -userName Webin-1 \
            -password password \
            -submit)

        # Check if the submission was successful
        if echo "$output" | grep -q "The submission has been completed successfully. The following run accession was assigned to the submission:"; then
            # Extract the run accession number
            run_accession=$(echo "$output" | grep "Run accession:" | awk '{print $NF}')
            
            # Log the successful submission with manifest and run accession
            echo "Submission successful for manifest: ${manifest##*/}, Run Accession: ${run_accession}" >> "$success_log"
        else
            # Log failed submissions
            echo "Submission FAILED for manifest: ${manifest##*/}" >> "$success_log"
            
            # Optionally, print the output for debugging
            echo "Submission output for ${manifest##*/}:" >> "$success_log"
            echo "$output" >> "$success_log"
            echo "---" >> "$success_log"
        fi
    done 
done