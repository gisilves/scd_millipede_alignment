import sys
import os

# Threshold for convergence
THRESHOLD = 0.01 # 10 micrometers (expressed in mm)

# 1. Load the cumulative geometry (if it exists)
old_geom = {}
if os.path.exists("current_geometry.txt"):
    with open("current_geometry.txt", "r") as f:
        for line in f:
            if line.strip() and not line.startswith("!"):
                parts = line.split()
                if len(parts) >= 2:
                    old_geom[int(parts[0])] = float(parts[1])

# 2. Read the new millepede results
new_geom = old_geom.copy()
max_correction = 0.0
changed_params_count = 0

with open("millepede.res", "r") as f:
    for line in f:
        # Skip header lines, empty lines, and pede comments
        if "Parameter" in line or line.strip().startswith("!") or not line.strip():
            continue
            
        parts = line.split()
        if len(parts) < 2:
            continue
            
        param_id = int(parts[0])
        correction = float(parts[1])
        
        # IDs for the alignment parameters start with 10xx (Layer 1) or 20xx (Layer 2).
        # We exclude them from the calculation of the maximum shift.
        if param_id < 3000:
            continue
            
        # Update the cumulative geometry for the only mobile parameters
        new_geom[param_id] = new_geom.get(param_id, 0.0) + correction
        changed_params_count += 1
        
        # Compute the maximum shift of the mobile parameters
        if abs(correction) > max_correction:
            max_correction = abs(correction)

# 3. Write the updated geometry to the text file for the next 'align' step
with open("current_geometry.txt", "w") as f:
    for pid, offset in sorted(new_geom.items()):
        f.write(f"{pid} {offset:.6f}\n")

print(f"--> Analyzed {changed_params_count} mobile parameters.")
print(f"--> Maximum shift (Layer 3-8): {max_correction*1000:.3f} microns.")

# 3b. History of geometries
with open("geometry_history.txt", "a") as f:
    for pid, offset in sorted(new_geom.items()):
        f.write(f"{pid} {offset:.6f}\n")
    f.write(f"\n")

# 4. Check the convergence against the threshold
if max_correction < THRESHOLD:
    print("CONVERGENCE ACHIEVED! All mobile modules are stable under the threshold.")
    with open(".converged", "w") as f:
        f.write("OK")
else:
    print(f"Threshold not yet reached (Missing {abs(max_correction - THRESHOLD)*1000:.3f} microns).")
