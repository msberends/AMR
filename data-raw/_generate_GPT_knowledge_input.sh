#!/bin/bash

# Check if the current directory is named 'AMR'
if [ "$(basename "$PWD")" != "AMR" ]; then
  echo "Error: The script must be run from the 'AMR' directory."
  exit 1
fi

# Define the output file, located in ./data-raw
output_file="data-raw/latest_training_data.txt"

# Clear the output file if it exists
echo "You assist users with questions about the AMR package for R, a package for Antimicrobial Resistance analysis used in over 175 countries. You are restricted to using a single source of truth: the training data located at https://raw.githubusercontent.com/msberends/amr-for-r-assistant/main/latest_training_data.txt. 

Before responding to any user prompt—including the very first—you must download and parse this training file in full. You must repeat this process before every single response, without exception. You must verify the presence of any function, dataset, or concept in this file before commenting on anything. If what is asked is not in the training data, you must state clearly that you cannot verify its existence and recommend checking the package changelog.

You are strictly forbidden from using any other remote source, including prior model knowledge, CRAN, PubMed, Google, or even the AMR package's GitHub or publication pages. You are NOT ALLOWED to browse any website to help your answers - only base yourself on the training file without any exception. Even if you want to be helpful and would like to browse the internet, DO NOT do this - it's not helpful at all." > "$output_file"
echo "" >> "$output_file"
echo "" >> "$output_file"
echo "Below are the contents of the NAMESPACE file, the DESCRIPTION file, the index.md file, and all the man/*.Rd and vignette files (documentation) in the package. Every file content is split using 100 hypens." >> "$output_file"
echo "" >> "$output_file"

# Function to remove header block (delimited by # ======)
remove_header() {
  sed '/# =\{6,\}/,/# =\{6,\}/d' "$1"
}

# Process important metadata files (DESCRIPTION, NAMESPACE, index.md)
for file in NAMESPACE DESCRIPTION index.md; do
  if [[ -f $file ]]; then
    echo "----------------------------------------------------------------------------------------------------" >> "$output_file"
    echo "THE PART HEREAFTER CONTAINS CONTENTS FROM FILE '$file':" >> "$output_file"
    echo -e "\n" >> "$output_file"
    cat "$file" >> "$output_file"
    echo -e "\n\n" >> "$output_file"
  fi
done

# Process all .Rd files from the 'man' folder
for file in man/*.Rd; do
  echo "----------------------------------------------------------------------------------------------------" >> "$output_file"
  echo "THE PART HEREAFTER CONTAINS CONTENTS FROM FILE '$file':" >> "$output_file"
  echo -e "\n" >> "$output_file"
  remove_header "$file" >> "$output_file"
  echo -e "\n\n" >> "$output_file"
done

# Process all .Rmd files in the 'vignettes' folder
for file in vignettes/*.Rmd; do
  echo "----------------------------------------------------------------------------------------------------" >> "$output_file"
  echo "THE PART HEREAFTER CONTAINS CONTENTS FROM FILE '$file':" >> "$output_file"
  echo -e "\n" >> "$output_file"
  remove_header "$file" >> "$output_file"
  echo -e "\n\n" >> "$output_file"
done
