#!/bin/bash

# Clean up old convergence files
rm -f .converged
rm -f current_geometry.txt

ITERATION=1
MAX_ITERATION=10

echo "=== STARTING ALIGNMENT LOOP ==="

while [ ! -f .converged ]
do
    echo ""
    echo "--------------------------------------------------"
    echo "Starting Iteration $ITERATION"
    echo "--------------------------------------------------"

    # 1. Run the local fit and generate the updated binary file
    echo "Running ./align..."
    ./align
    
    # 2. Start Millepede II to calculate global corrections
    echo "Running pede..."
    ./target/pede steer.txt > /dev/null
    
    if [ ! -f millepede.res ]; then
        echo "Error: millepede.res not generated. Exiting loop."
        exit 1
    fi

    # 3. Process the results and check the convergence threshold
    python3 update_geometry.py

    # Safety check to avoid infinite loops
    if [ $ITERATION -ge $MAX_ITERATION ]; then
        echo "Max iteration reached without convergence."
        break
    fi

    ((ITERATION++))
done

echo "=== END OF ALIGNMENT LOOP ==="