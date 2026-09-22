import pandas as pd
import numpy as np
import argparse


parser = argparse.ArgumentParser(description="Apply artificial misalignments to a hits CSV file.")
parser.add_argument("input", help="Input hits CSV file")
parser.add_argument("output", help="Output hits CSV file")
args = parser.parse_args()

# ==============================================================================
# 1. CONFIGURATION OF ARTIFICIAL MISALIGNMENTS (Known values for testing)
# ==============================================================================
# Nominal Z positions of layers (in mm)
Z_LAYERS = [0.0, 5.0, 55.0, 60.0, 110.0, 115.0, 165.0, 170.0]
IS_X_LAYER = [False, True, False, True, False, True, False, True]
IS_ACTIVE = [True, False, True, True, True, True, False, True]

# Define the coupled displacements to inject (in mm or radians)
MISALIGNMENTS = {
    # Layer: [Shift_X/Y, Rot_Z, Tilt_X, Tilt_Y]
    0: [+0.00,     0.0,     0.0,    0.0],     # Layer 0 (Y, REFERENCE - fixed)
    1: [+0.00,     0.0,     0.0,    0.0],     # Layer 1 (X, inactive)
    2: [+0.23,     0.0,     0.0,    0.0],     # Layer 2 (Y, free)
    3: [+0.00,     0.0,     0.0,    0.0],     # Layer 3 (X, REFERENCE - fixed)
    4: [-0.40,     0.0,     0.0,    0.0],     # Layer 4 (Y, free)
    5: [+0.88,     0.0,     0.0,    0.0],     # Layer 5 (X, free)
    6: [+0.00,     0.0,     0.0,    0.0],     # Layer 6 (Y, inactive)
    7: [+0.14,     0.0,     0.0,    0.0],     # Layer 7 (X, free)
}

# ==============================================================================
# 2. GEOMETRIC TRANSFORMATION PROCESS
# ==============================================================================
file_input = args.input 
file_output = args.output

# Load original hits
df = pd.read_csv(file_input, header=None)
df_misaligned = df.copy()

print(f"Loaded {len(df)} tracks from {file_input}")
print("Applying artificial misalignments...")

# For each track, extract an estimate of the trajectory (x,y) at that Z position
# to correctly apply rotations and tilts, which depend on impact coordinates.
for idx, row in df.iterrows():
    hits = row.to_numpy()
    
    # On-the-fly fit of ideal track to know x(z) and y(z) for this specific particle
    z_x = [Z_LAYERS[i] for i in range(8) if IS_X_LAYER[i]]
    h_x = [hits[i] for i in range(8) if IS_X_LAYER[i]]
    z_y = [Z_LAYERS[i] for i in range(8) if not IS_X_LAYER[i]]
    h_y = [hits[i] for i in range(8) if not IS_X_LAYER[i]]
    
    tx, x0 = np.polyfit(z_x, h_x, 1)
    ty, y0 = np.polyfit(z_y, h_y, 1)
    
    # Apply deformations layer by layer
    for i in range(8):
        z = Z_LAYERS[i]
        x_ideal = x0 + tx * z
        y_ideal = y0 + ty * z
        
        # Extract injected deltas for this specific layer
        shift, rot_z, tilt_x, tilt_y = MISALIGNMENTS[i]
        
        if IS_X_LAYER[i]:
            # Geometric effect of deformations on X coordinate read by sensor:
            # delta_x = Shift_X - y*Rot_Z + y*tx*Tilt_X - z*tx*Tilt_Y
            delta = shift - (y_ideal * rot_z) + (y_ideal * tx * tilt_x) - (z * tx * tilt_y)
            df_misaligned.iat[idx, i] += delta
        else:
            # Geometric effect of deformations on Y coordinate read by sensor:
            # delta_y = Shift_Y + x*Rot_Z + z*ty*Tilt_X - x*ty*Tilt_Y
            delta = shift + (x_ideal * rot_z) + (z * ty * tilt_x) - (x_ideal * ty * tilt_y)
            df_misaligned.iat[idx, i] += delta

# Save new altered CSV
df_misaligned.to_csv(file_output, header=False, index=False)

print(f"Done. Misaligned file saved to: {file_output}")
print("\nExpected values that Millepede SHOULD find:")
for k, v in MISALIGNMENTS.items():
    if any(val != 0.0 for val in v):
        print(f"  Layer {k+1} -> Shift: {v[0]:.3f} mm, RotZ: {v[1]:.4f} rad")