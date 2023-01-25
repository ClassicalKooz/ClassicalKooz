#!/bin/bash

#Loop through 10 subdirectories
for i in $(seq 1 10); do
  #Change into the i-th subdirectory
  cd ./replicate$i

  #Find all files with "model4" in their name
  files=$(find . -name "*model4*")

  #Loop through all matching files
  for file in $files; do
    #Rename the file from model4 to model2
    mv "$file" "${file/model7/model3}"
  done

  #Change back to the parent directory
  cd ..
done