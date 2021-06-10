
import os
import numpy as np
from pygeosx_tools import wrapper, parallel_io, plot_tools
import hdf5_wrapper
import matplotlib.pyplot as plt
from matplotlib import cm
from pyevtk.hl import gridToVTK
from scipy import interpolate


class InSAR_Analysis():
    def __init__(self):
        """
        InSAR Analysis class
        """
        self.set_names = []
        self.set_keys = []
        self.local_insar_map = []

        self.node_position_key = ''
        self.node_displacement_key = ''
        self.node_ghost_key = ''

        self.satellite_vector = [0.0, 0.0, 1.0]
        self.wavelength = 1.0

        self.x_grid = []
        self.y_grid = []
        self.elevation = 0.0
        self.times = []
        self.mask = []
        self.displacement = []
        self.range_change = []
        self.phase = []

    def setup_grid(self, problem, set_names=[], x_range=[], y_range=[], dx=1.0, dy=1.0):
        """
        Setup the InSAR grid

        Args:
            problem (pygeosx.group): GEOSX problem handle
            set_names (list): List of set names to apply to the analysis to
            x_range (list): Extents of the InSAR image in the x-direction (optional)
            y_range (list): Extents of the InSAR image in the y-direction (optional)
            dx (float): Resolution of the InSAR image in the x-direction (default=1)
            dy (float): Resolution of the InSAR image in the y-direction (default=1)
        """
        # Determine pygeosx keys
        self.set_names = set_names
        self.set_keys = [wrapper.get_matching_wrapper_path(problem, ['nodeManager', s]) for s in set_names]
        self.node_position_key = wrapper.get_matching_wrapper_path(problem, ['nodeManager', 'ReferencePosition'])
        self.node_displacement_key = wrapper.get_matching_wrapper_path(problem, ['nodeManager', 'TotalDisplacement'])
        self.node_ghost_key = wrapper.get_matching_wrapper_path(problem, ['nodeManager', 'ghostRank'])

        # If not specified, then setup the grid extents
        ghost_rank = wrapper.get_wrapper(problem, self.node_ghost_key)
        x = wrapper.get_wrapper(problem, self.node_position_key)
        set_ids = np.concatenate([wrapper.get_wrapper(problem, k) for k in self.set_keys], axis=0)

        # Choose non-ghost set members
        xb = x[set_ids, :]
        gb = ghost_rank[set_ids]
        xc = xb[gb < 0, :]
        global_min, global_max = parallel_io.get_global_array_range(xc)

        if (len(x_range) == 0):
            x_range = [global_min[0], global_max[0]]
        if (len(y_range) == 0):
            y_range = [global_min[1], global_max[1]]

        # Choose the grid
        Nx = int(np.ceil((x_range[1] - x_range[0]) / dx))
        Ny = int(np.ceil((y_range[1] - y_range[0]) / dy))
        self.x_grid = np.linspace(x_range[0], x_range[1], Nx + 1)
        self.y_grid = np.linspace(y_range[0], y_range[1], Ny + 1)

        # Save the average elevation for vtk outputs
        self.elevation = global_min[0]

        # Trigger the map build
        self.build_map(problem)

    def build_map(self, problem):
        """
        Build the map between the mesh, InSAR image.
        Note: this method can be used to update the map after
              significant changes to the mesh.

        Args:
            problem (pygeosx.group): GEOSX problem handle
        """
        # Load the data
        ghost_rank = wrapper.get_wrapper(problem, self.node_ghost_key)
        x = wrapper.get_wrapper(problem, self.node_position_key)
        set_ids = np.concatenate([wrapper.get_wrapper(problem, k) for k in self.set_keys], axis=0)

        # Choose non-ghost set members
        xb = x[set_ids, :]
        gb = ghost_rank[set_ids]
        xc = xb[gb < 0, :]

        # Setup the node to insar map
        self.local_insar_map = []
        dx = self.x_grid[1] - self.x_grid[0]
        dy = self.y_grid[1] - self.y_grid[0]
        x_bins = np.concatenate([[self.x_grid[0] - 0.5 * dx],
                                 0.5*(self.x_grid[1:] + self.x_grid[:-1]),
                                 [self.x_grid[-1] + 0.5 * dx]])
        y_bins = np.concatenate([[self.y_grid[0] - 0.5 * dy],
                                 0.5*(self.y_grid[1:] + self.y_grid[:-1]),
                                 [self.y_grid[-1] + 0.5 * dy]])
        if len(xc):
            Ix = np.digitize(np.squeeze(xc[:, 0]), x_bins) - 1
            Iy = np.digitize(np.squeeze(xc[:, 1]), y_bins) - 1
            for ii in range(len(self.x_grid)):
                for jj in range(len(self.y_grid)):
                    tmp = np.where((Ix == ii) & (Iy == jj))[0]
                    if len(tmp):
                        self.local_insar_map.append([ii, jj, tmp])

    def extract_insar(self, problem):
        """
        Extract InSAR image for current step

        Args:
            problem (pygeosx.group): GEOSX problem handle
        """
        # Load values
        time = wrapper.get_wrapper(problem, 'Events/time')[0]
        ghost_rank = wrapper.get_wrapper(problem, self.node_ghost_key)
        x = wrapper.get_wrapper(problem, self.node_displacement_key)
        set_ids = np.concatenate([wrapper.get_wrapper(problem, k) for k in self.set_keys], axis=0)

        # Choose non-ghost set members
        xb = x[set_ids, :]
        gb = ghost_rank[set_ids]
        xc = xb[gb < 0, :]

        # Find local displacements
        Nx = len(self.x_grid)
        Ny = len(self.y_grid)
        local_displacement_sum = np.zeros((Nx, Ny, 3))
        local_N = np.zeros((Nx, Ny), dtype='int')
        for m in self.local_insar_map:
            local_N[m[0], m[1]] += len(m[2])
            for ii in m[2]:
                local_displacement_sum[m[0], m[1], :] += xc[ii, :]

        # Communicate values
        global_displacement_sum = np.sum(np.array(parallel_io.gather_array(local_displacement_sum, concatenate=False)), axis=0)
        global_N = np.sum(np.array(parallel_io.gather_array(local_N, concatenate=False)), axis=0)

        # Find final 3D displacement
        global_displacement = np.zeros((Nx, Ny, 3))
        if (parallel_io.rank == 0):
            range_change = np.zeros((Nx, Ny))
            for ii in range(3):
                d = np.squeeze(global_displacement_sum[:, :, ii]) / (global_N + 1e-10)
                d[global_N == 0] = np.NaN
                global_displacement[:, :, ii] = self.fill_nan_gaps(d)
                range_change += d * self.satellite_vector[ii]

            # Filter nans
            self.displacement.append(global_displacement)
            self.range_change.append(range_change)
            self.phase.append(np.angle(np.exp(2j * np.pi * range_change / self.wavelength)))
            self.mask.append(global_N > 0)
            self.times.append(time)

    def fill_nan_gaps(self, values):
        """
        Fill gaps in the insar data which are specified via NaN's
        """
        z = np.isnan(values)
        if np.sum(z):
            N = np.shape(values)
            grid = np.meshgrid(self.x_grid, self.y_grid, indexing='ij')

            # Filter out values in flattened arrays
            x_flat = np.reshape(grid[0], (-1))
            y_flat = np.reshape(grid[1], (-1))
            v_flat = np.reshape(values, (-1))
            t_flat = np.reshape(z, (-1))
            I_valid = np.where(~t_flat)[0]
            x_valid = x_flat[I_valid]
            y_valid = y_flat[I_valid]
            v_valid = v_flat[I_valid]

            # Re-interpolate values
            vb = interpolate.griddata((x_valid, y_valid), v_valid, tuple(grid), method='linear')
            values = np.reshape(vb, N)
        return values

    def save_hdf5(self, header='insar', output_root='./results'):
        os.makedirs(output_root, exist_ok=True)
        with hdf5_wrapper.hdf5_wrapper('%s/%s.hdf5' % (output_root, header), mode='w') as data:
            data['x'] = self.x_grid
            data['y'] = self.y_grid
            data['time'] = self.times
            data['displacement'] = np.array(self.displacement)
            data['range_change'] = np.array(self.range_change)
            data['phase'] = np.array(self.phase)

    def save_csv(self, header='insar', output_root='./results'):
        os.makedirs(output_root, exist_ok=True)
        np.savetxt('%s/%s_x_grid.csv' % (output_root, header),
                   self.x_grid,
                   delimiter=', ')
        np.savetxt('%s/%s_y_grid.csv' % (output_root, header),
                   self.y_grid,
                   delimiter=', ')
        for ii, t in enumerate(self.times):
            comments = 'T (days), %1.8e' % (t / (60 * 60 * 24))
            np.savetxt('%s/%s_range_change_%03d.csv' % (output_root, header, ii),
                       self.range_change[ii],
                       delimiter=', ',
                       header=comments)
            np.savetxt('%s/%s_phase_%03d.csv' % (output_root, header, ii),
                       self.phase[ii],
                       delimiter=', ',
                       header=comments)

    def save_vtk(self, header='insar', output_root='./results'):
        os.makedirs(output_root, exist_ok=True)
        x = np.ascontiguousarray(self.x_grid)
        y = np.ascontiguousarray(self.y_grid)
        z = np.array([self.elevation])

        for ii, t in enumerate(self.times):
            data = {'range_change': np.ascontiguousarray(np.expand_dims(self.range_change[ii], -1)),
                    'phase': np.ascontiguousarray(np.expand_dims(self.phase[ii], -1)),
                    'dx': np.ascontiguousarray(np.expand_dims(self.displacement[ii][:, :, 0], -1)),
                    'dy': np.ascontiguousarray(np.expand_dims(self.displacement[ii][:, :, 1], -1)),
                    'dz': np.ascontiguousarray(np.expand_dims(self.displacement[ii][:, :, 2], -1))}

            gridToVTK('%s/%s_%03d' % (output_root, header, ii),
                      x,
                      y,
                      z,
                      pointData=data)

    def save_image(self, header='insar', output_root='./results', interp_method='quadric'):
        os.makedirs(output_root, exist_ok=True)
        fig = plot_tools.HighResPlot()

        for ii, t in enumerate(self.times):
            # Range change
            fig.reset()
            extents = [self.x_grid[0], self.x_grid[-1], self.y_grid[0], self.y_grid[-1]]
            ca = plt.imshow(np.transpose(np.flipud(self.range_change[ii])),
                            extent=extents,
                            cmap=cm.jet,
                            aspect='auto',
                            interpolation=interp_method)
            plt.title('T = %1.4e (days)' % (t / (60 * 60 * 24)))
            plt.xlabel('X (m)')
            plt.ylabel('Y (m)')
            cb = plt.colorbar(ca)
            cb.set_label('Range Change (m)')
            fig.save('%s/%s_range_change_%03d' % (output_root, header, ii))

            # Wrapped phase
            fig.reset()
            extents = [self.x_grid[0], self.x_grid[-1], self.y_grid[0], self.y_grid[-1]]
            ca = plt.imshow(np.transpose(np.flipud(self.phase[ii])),
                            extent=extents,
                            cmap=cm.jet,
                            aspect='auto',
                            interpolation=interp_method)
            plt.title('T = %1.4e (days)' % (t / (60 * 60 * 24)))
            plt.xlabel('X (m)')
            plt.ylabel('Y (m)')
            cb = plt.colorbar(ca)
            cb.set_label('Phase (radians)')
            fig.save('%s/%s_wrapped_phase_%03d' % (output_root, header, ii))


