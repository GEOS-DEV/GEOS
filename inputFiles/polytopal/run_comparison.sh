
#!/bin/bash

# It's set to exit immediately if any command fails (set -e).
set -e

echo "Starting GEOSX batch simulations..."

# --- First Simulation ---
echo "Running simulation with tpfa.xml input..."
../../build-macOS_arm-release/bin/Debug/geosx -i buoyancy_column_co2_3d_tpfa.xml -o ../../../geos_output
echo "First simulation finished."

# --- Second Simulation ---
echo "Running simulation with mfd.xml input..."
../../build-macOS_arm-release/bin/Debug/geosx -i buoyancy_column_co2_3d_mfd.xml -o ../../../geos_output
echo "Second simulation finished."

echo "All simulations completed successfully."

