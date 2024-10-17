#!/bin/bash

# Check if the current directory is named 'AMR'
if [ "$(basename "$PWD")" != "AMR" ]; then
  echo "Error: The script must be run from the 'AMR' directory."
  exit 1
fi

# Define the output file, located in ./data-raw
output_file="data-raw/gpt_training_text.txt"

# Clear the output file if it exists
echo "This files contains all context you must know about the AMR package for R."> "$output_file"
echo -e "\n\n\n\n" >> "$output_file"

# Function to remove header block (delimited by # ======)
remove_header() {
  sed '/# =\{6,\}/,/# =\{6,\}/d' "$1"
}

# Process all .R files in the 'R' folder
for file in R/*.R; do
  echo "THE NEXT PART CONTAINS CONTENTS FROM FILE $file" >> "$output_file"
  echo -e "\n\n" >> "$output_file"
  remove_header "$file" >> "$output_file"
  echo -e "\n\n" >> "$output_file"
done

# Process all .Rmd files in the 'vignettes' folder
for file in vignettes/*.Rmd; do
  echo "THE NEXT PART CONTAINS CONTENTS FROM FILE $file" >> "$output_file"
  echo -e "\n\n" >> "$output_file"
  remove_header "$file" >> "$output_file"
  echo -e "\n\n" >> "$output_file"
done

# Process important metadata files (DESCRIPTION, NAMESPACE, index.md)
for file in DESCRIPTION NAMESPACE index.md; do
  if [[ -f $file ]]; then
    echo "THE NEXT PART CONTAINS CONTENTS FROM FILE $file" >> "$output_file"
    echo -e "\n\n" >> "$output_file"
    cat "$file" >> "$output_file"
    echo -e "\n\n" >> "$output_file"
  fi
done

# Process test files (if available) in the 'tests' folder
for file in tests/*.R; do
  echo "THE NEXT PART CONTAINS CONTENTS FROM FILE $file" >> "$output_file"
  echo -e "\n\n" >> "$output_file"
  remove_header "$file" >> "$output_file"
  echo -e "\n\n" >> "$output_file"
done

# Process all .Rd files from the 'man' folder
for file in man/*.Rd; do
  echo "THE NEXT PART CONTAINS CONTENTS FROM FILE $file" >> "$output_file"
  echo -e "\n\n" >> "$output_file"
  remove_header "$file" >> "$output_file"
  echo -e "\n\n" >> "$output_file"
done
