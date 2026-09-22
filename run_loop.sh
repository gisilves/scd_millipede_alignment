#!/bin/bash

# Clean up old convergence files
rm -f .converged
rm -f current_geometry.txt
rm -f rootfiles/*
rm -f millepede.*
rm -f mille_output_run001.bin
rm -f geometry_history.txt

ITERATION=1
MAX_ITERATION=10

echo "=== STARTING ALIGNMENT LOOP ==="

# Check if the first argument is provided
if [ -z "$1" ]; then
  echo "No argument provided. Using default file "hits_misaligned.csv""
  FILE="hits_misaligned.csv"
else
  FILE=$1
fi

while [ ! -f .converged ]
do
    echo ""
    echo "--------------------------------------------------"
    echo "Starting Iteration $ITERATION"
    echo "--------------------------------------------------"

    # 1. Run the local fit and generate the updated binary file
    echo "Running ./align..."
    ./align $ITERATION $FILE
    
    # 2. Start Millepede II to calculate global corrections
    echo "Running pede..."
    ./target/pede steer.txt >> /dev/null
    
    if [ ! -f millepede.res ]; then
        echo "Error: millepede.res not generated. Exiting loop."
        exit 1
    fi

    # 3. Process the results and check the convergence threshold
    python3 python/update_geometry.py

    # Safety check to avoid infinite loops
    if [ $ITERATION -ge $MAX_ITERATION ]; then
        echo "Max iteration reached without convergence."
        break
    fi

    ((ITERATION++))
done

echo "=== END OF ALIGNMENT LOOP ==="