#!/usr/bin/env bash
set -euo pipefail
python3 generateRigidBodyJamming2DParticles.py
geosx -i mpm_rigidBodyJamming2D.xml
python3 checkRigidBodyJamming2D.py
