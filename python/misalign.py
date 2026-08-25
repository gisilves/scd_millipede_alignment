import pandas as pd
import numpy as np

# ==============================================================================
# 1. CONFIGURATION OF ARTIFICIAL MISALIGNMENTS (Known values for testing)
# ==============================================================================
# Nominal Z positions of layers (in mm)
Z_LAYERS = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0]
IS_X_LAYER = [True, False, True, False, True, False, True, False]

# Define the coupled displacements to inject (in mm or radians)
MISALIGNMENTS = {
    # Layer: [Shift_X/Y, Rot_Z, Tilt_X, Tilt_Y]
    0: [0.0,     0.0,     0.0,    0.0],     # Layer 1 (X)
    1: [0.0,     0.0,     0.0,    0.0],     # Layer 2 (Y)
    2: [0.0,     0.0,     0.0,    0.0],     # Layer 3 (X)
    3: [0.0,     0.0,     0.0,    0.0],     # Layer 4 (Y)
    4: [0.0,     0.0,     0.0,    0.0],     # Layer 5 (X)
    5: [0.0,     0.0,     0.0,    0.0],     # Layer 6 (Y)
    6: [0.0,     0.0,     0.0,    0.0],     # Layer 7 (X)
    7: [0.0,     0.0,     0.0,    0.0],     # Layer 8 (Y)
}

# ==============================================================================
# 2. GEOMETRIC TRANSFORMATION PROCESS
# ==============================================================================
file_input = "hits_perfect.csv"
file_output = "hits_generated.csv"

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