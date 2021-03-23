#ifndef GRID_H
#define GRID_H

#include "constants.hpp"

struct grid_t {
    llint ntaperx, ntapery, ntaperz;
    llint ndampx, ndampy, ndampz;
    llint nx, ny, nz;
    llint dx, dy, dz;
    llint x1, x2, x3, x4, x5, x6;
    llint y1, y2, y3, y4, y5, y6;
    llint z1, z2, z3, z4, z5, z6;
    llint lx, ly, lz;
    // These parameters are used only for tasks
    llint ntx, nty, tsx, tsy;

    // Used for multi-gpu targets
    llint ngpu;
};

struct grid_t init_grid(llint nx, llint ny, llint nz, llint tsx, llint tsy,
                        llint ngpu);

#endif
