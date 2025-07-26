# Active Learning

Scripts for analyzing uncertainty in MD trajectories and selecting data for retraining.

## Contents

- Identify atom-step combinations from MD data (stored in HDF5) where local uncertainty exceeds a threshold ("spikes").
  - Outputs: all atom-steps, unique spike steps, total count, and first 100 unique steps.

- Compute average spike thresholds across bonding regimes:
  - Surface Pt, Intermediate Pt, Surface H, Gas phase H.
  - Select 100 representative spike-containing frames for retraining.

