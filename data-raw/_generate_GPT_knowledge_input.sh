#!/bin/bash

# Check if the current directory is named 'AMR'
if [ "$(basename "$PWD")" != "AMR" ]; then
  echo "Error: The script must be run from the 'AMR' directory."
  exit 1
fi

rm -rf data-raw/gpt_training_text_v*

# Define the output file, located in ./data-raw
version="$1"
output_file="data-raw/gpt_training_text_v${version}.txt"

# Clear the output file if it exists
echo "This knowledge base contains all context you must know about the AMR package for R. You are a GPT trained to be an assistant for the AMR package in R. You are an incredible R specialist, especially trained in this package and in the tidyverse." > "$output_file"
echo "" >> "$output_file"
echo "First and foremost, you are trained on version ${version}. Remember this whenever someone asks which AMR package version you’re at." >> "$output_file"
echo "" >> "$output_file"
echo "Below are the contents of the `NAMESPACE` file, the `index.md` file, and all the `man/*.Rd` files (documentation) in the package. Every file content is split using 100 hypens." >> "$output_file"
echo "----------------------------------------------------------------------------------------------------" >> "$output_file"
echo "" >> "$output_file"

# Function to remove header block (delimited by # ======)
remove_header() {
  sed '/# =\{6,\}/,/# =\{6,\}/d' "$1"
}

# # Process all .R files in the 'R' folder
# for file in R/*.R; do
#   echo "--------------------------------------------------" >> "$output_file"
#   echo "THE PART HEREAFTER CONTAINS CONTENTS FROM FILE '$file':" >> "$output_file"
#   echo -e "\n" >> "$output_file"
#   remove_header "$file" >> "$output_file"
#   echo -e "\n\n" >> "$output_file"
# done

# Process important metadata files (DESCRIPTION, NAMESPACE, index.md)
for file in NAMESPACE index.md; do
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

# Process README.md
# echo "THE PART HEREAFTER CONTAINS THE README OF OUR PYTHON PACKAGE" >> "$output_file"
# echo -e "\n" >> "$output_file"
# for file in PythonPackage/AMR/README.md; do
#   remove_header "$file" >> "$output_file"
#   echo -e "\n\n" >> "$output_file"
# done

# Process test files (if available) in the 'tests' folder
# for file in tests/*.R; do
#   echo "THE PART HEREAFTER CONTAINS CONTENTS FROM FILE '$file':" >> "$output_file"
#   echo -e "\n" >> "$output_file"
#   remove_header "$file" >> "$output_file"
#   echo -e "\n\n" >> "$output_file"
# done
