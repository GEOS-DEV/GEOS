#!/usr/bin/env python3
"""
Compare computed gravity data from HDF5 file to analytic rectangular prism solution.

Analytic solution reference:
Nagy, D. (1966). "The gravitational attraction of a right rectangular prism." Geophysics, 31(2), 362-371. https://doi.org/10.1190/1.1439779
"""

import h5py
import numpy as np
import sys

# Analytic solution for gravity anomaly (gz) due to a rectangular prism
def analyticPrismGz(station, prism_bounds, density_contrast):
    x, y, z = station
    x0, x1, y0, y1, z0, z1 = prism_bounds
    epsilon = 1e-16
    val = 0.0
    for i in range(2):
        for j in range(2):
            for k in range(2):
                xi = x0 if i == 0 else x1
                yj = y0 if j == 0 else y1
                zk = z0 if k == 0 else z1
                dx = xi - x
                dy = yj - y
                dz = zk - z
                r = np.sqrt(dx*dx + dy*dy + dz*dz)
                if r > epsilon:
                    sign = 1.0 if (i + j + k) % 2 == 0 else -1.0
                    arg1 = np.arctan2(dx * dy, dz * r + epsilon)
                    arg2 = np.log(r + dy + epsilon)
                    arg3 = np.log(r + dx + epsilon)
                    val += sign * (dz * arg1 - dx * arg2 - dy * arg3)
    G = 6.67430e-11  # m^3 kg^-1 s^-2
    return G * density_contrast * val


def write_vtk_surface(filename, x_unique, y_unique, z_val, gravity_2d, analytic_2d=None):
    """Write gravity data as a VTK RECTILINEAR_GRID surface."""
    nx, ny = len(x_unique), len(y_unique)
    with open(filename, 'w') as f:
        f.write('# vtk DataFile Version 3.0\n')
        f.write('Gravity surface data\n')
        f.write('ASCII\n')
        f.write(f'DATASET RECTILINEAR_GRID\n')
        f.write(f'DIMENSIONS {nx} {ny} 1\n')
        f.write(f'X_COORDINATES {nx} float\n')
        for xi in x_unique:
            f.write(f'{xi} ')
        f.write('\n')
        f.write(f'Y_COORDINATES {ny} float\n')
        for yi in y_unique:
            f.write(f'{yi} ')
        f.write('\n')
        f.write(f'Z_COORDINATES 1 float\n')
        f.write(f'{z_val}\n')
        f.write(f'\nPOINT_DATA {nx*ny}\n')
        f.write('SCALARS gravity float 1\n')
        f.write('LOOKUP_TABLE default\n')
        for iy in range(ny):
            for ix in range(nx):
                f.write(f'{gravity_2d[ix, iy]}\n')
        if analytic_2d is not None:
            f.write('\nSCALARS analytic float 1\n')
            f.write('LOOKUP_TABLE default\n')
            for iy in range(ny):
                for ix in range(nx):
                    f.write(f'{analytic_2d[ix, iy]}\n')


def compare_gravity_data(hdf5_file, prism_bounds, density_contrast, save_fig=None):
    """Extract gravity data and coordinates from HDF5 file and compare to analytic solution"""
    try:
        with h5py.File(hdf5_file, 'r') as f:
            # Read coordinate data
            coord_x = f['/coordinateX'][:]
            coord_y = f['/coordinateY'][:]
            coord_z = f['/coordinateZ'][:]
            # Read gravity data
            gravity_data = f['/scatterData'][:]
            # Read time data
            time_data = f['/gravity Time'][:]
            # Assume single time step for simplicity
            if len(coord_x.shape) > 1:
                current_x = coord_x[0, :]
                current_y = coord_y[0, :]
                current_z = coord_z[0, :]
                current_gravity = gravity_data[0, :]
            else:
                current_x = coord_x
                current_y = coord_y
                current_z = coord_z
                current_gravity = gravity_data
            stations = np.stack([current_x, current_y, -current_z], axis=-1)
            # Compute analytic solution
            gz_analytic = np.array([analyticPrismGz(station, prism_bounds, density_contrast) for station in stations])
            # Compare
            diff = current_gravity - gz_analytic
            for i, (g_comp, g_ana, d) in enumerate(zip(current_gravity, gz_analytic, diff)):
                print(f"Station {i:4d}: Computed={g_comp:15.6e} Analytic={g_ana:15.6e} Diff={d:15.6e}")
            print("Comparison to analytic solution:")
            print(f"  Max abs difference: {np.max(np.abs(diff)):.6e}")
            print(f"  Mean abs difference: {np.mean(np.abs(diff)):.6e}")

            # Plot results as 2D arrays
            try:
                import matplotlib.pyplot as plt
                # Infer grid shape from coordinate data
                unique_x = np.sort(np.unique(current_x))
                unique_y = np.sort(np.unique(current_y))
                nx = unique_x.size
                ny = unique_y.size
                # Build 2D arrays using coordinate mapping
                gz_analytic_2d = np.zeros((nx, ny))
                current_gravity_2d = np.zeros((nx, ny))
                diff_2d = np.zeros((nx, ny))
                x_grid = np.zeros((nx, ny))
                y_grid = np.zeros((nx, ny))
                # Map each station to grid index
                for idx in range(len(current_x)):
                    ix = np.searchsorted(unique_x, current_x[idx])
                    iy = np.searchsorted(unique_y, current_y[idx])
                    gz_analytic_2d[ix, iy] = gz_analytic[idx]
                    current_gravity_2d[ix, iy] = current_gravity[idx]
                    diff_2d[ix, iy] = diff[idx]
                    x_grid[ix, iy] = current_x[idx]
                    y_grid[ix, iy] = current_y[idx]
                rel_diff_2d = 100.0 * np.abs(diff_2d / (gz_analytic_2d + 1e-16))  # percentage error
                # Use min/max for extent
                extent = [y_grid.min(), y_grid.max(), x_grid.min(), x_grid.max()]

                # Write VTK surface output
                z_val = current_z[0] if len(current_z) > 1 else current_z
                write_vtk_surface('gravity_surface.vtk', unique_x, unique_y, z_val, current_gravity_2d, analytic_2d=gz_analytic_2d)

                # Set common colorbar range for analytic and computed
                vmin = min(gz_analytic_2d.min(), current_gravity_2d.min())
                vmax = max(gz_analytic_2d.max(), current_gravity_2d.max())
                fig, axs = plt.subplots(1, 3, figsize=(18, 6))

                # Prism outline coordinates (define here before plotting)
                x0, x1, y0, y1, z0, z1 = prism_bounds
                prism_x = [y0, y1, y1, y0, y0]
                prism_y = [x0, x0, x1, x1, x0]

                im0 = axs[0].imshow(gz_analytic_2d, aspect='equal', cmap='viridis', extent=extent, origin='lower', vmin=vmin, vmax=vmax)
                axs[0].set_title('Analytic Solution')
                plt.colorbar(im0, ax=axs[0],  label='gz (m/s²)')
                axs[0].plot(prism_x, prism_y, color='red', linewidth=2, label='Prism Outline')

                im1 = axs[1].imshow(current_gravity_2d, aspect='equal', cmap='viridis', extent=extent, origin='lower', vmin=vmin, vmax=vmax)
                axs[1].set_title('Computed (HDF5)')
                plt.colorbar(im1, ax=axs[1], label='gz (m/s²)')
                axs[1].plot(prism_x, prism_y, color='red', linewidth=2, label='Prism Outline')

                im2 = axs[2].imshow(rel_diff_2d, aspect='equal', cmap='magma', extent=extent, origin='lower')
                axs[2].set_title('Relative Difference (%)')
                cbar = plt.colorbar(im2, ax=axs[2])
                cbar.set_label('Percent Error (%)')
                axs[2].plot(prism_x, prism_y, color='red', linewidth=2, label='Prism Outline')

                for ax in axs:
                    ax.set_xlabel('Y Coordinate')
                    ax.set_ylabel('X Coordinate')
                    ax.set_aspect('equal')
                    ax.legend()
                plt.suptitle('Gravity Comparison: Analytic vs Computed')
                plt.tight_layout(rect=[0, 0.03, 1, 0.95])
                if save_fig:
                    plt.savefig(save_fig, dpi=300)
                    print(f"Figure saved to {save_fig}")
                plt.show()
                
                # Numerical accuracy test (after plotting)
                max_error = np.max(rel_diff_2d)
                mean_error = np.mean(rel_diff_2d)
                threshold_max = 1.0  # percent
                print(f"Max percent error: {max_error:.4f}%")
                passed = (max_error < threshold_max)
                if passed:
                    print("Numerical accuracy test: PASS")
                else:
                    print("Numerical accuracy test: FAIL")
                return passed
            except ImportError:
                print("matplotlib not installed. Skipping plot.")
                return False
    except Exception as e:
        print(f"Error reading HDF5 file: {e}")
        return False

if __name__ == "__main__":
    if len(sys.argv) < 2:
        input_file = "gravity_history.hdf5"
    else:
        input_file = sys.argv[1]
    density_contrast = 900.0
    prism_bounds = [10000, 14000, 8000, 11000, 2000, 2200]
    
    print(f"Comparing gravity data from: {input_file}")
    print(f"Density contrast: {density_contrast}")
    print(f"Prism bounds: {prism_bounds}")
    save_fig = "gravity_comparison.png"
    success = compare_gravity_data(input_file, prism_bounds, density_contrast, save_fig)
    if success:
        print("Analysis completed successfully.")
    else:
        print("Analysis failed.")
