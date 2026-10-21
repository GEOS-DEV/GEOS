# Native geometry has moved to the unified representation

See [UNIFIED_PARTICLE_GEOMETRY.md](UNIFIED_PARTICLE_GEOMETRY.md) for the current
schema, minimal storage, reset policies, CPDI integration, CZ behavior, and tests.

The former optional glyph-reference matrix and legacy particle-file columns
have been removed. Regenerate all particle files and XML with the new PFW.
Old Python inputs using `particleReferenceGeometry` remain accepted as aliases.
Old GEOS input/restart files are intentionally unsupported.
