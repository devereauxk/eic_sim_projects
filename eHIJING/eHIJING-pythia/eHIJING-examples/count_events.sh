#!/bin/bash

# Set the directory path containing .dat files
directory=$1

# Initialize the sum
sum=0

# Loop through each .dat file in the directory
for file in "$directory"/*.dat; do
    # Extract the number from the last line before the first comma
    number=$(tail -n 1 "$file" | cut -d ',' -f 1)
    
    # Add the number to the sum
    sum=$((sum + number))
done

# Print the total sum
echo "Total events in $directory: $sum"
