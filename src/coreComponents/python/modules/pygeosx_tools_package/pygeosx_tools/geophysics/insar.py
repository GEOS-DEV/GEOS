
import os
import numpy as np
from pygeosx_tools import wrapper, parallel_io, plot_tools
import matplotlib.pyplot as plt
from matplotlib import cm


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
        self.satellite_wavelength = 1.0

        self.x_grid = []
        self.y_grid = []
        self.times = []
        self.displacement = []
        self.range_change = []

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

        # Setup the node to insar map
        x_bins = np.concatenate([[self.x_grid[0] - 0.5 * dx],
                                 0.5*(self.x_grid[1:] + self.x_grid[:-1]),
                                 [self.x_grid[-1] + 0.5 * dx]])
        y_bins = np.concatenate([[self.y_grid[0] - 0.5 * dy],
                                 0.5*(self.y_grid[1:] + self.y_grid[:-1]),
                                 [self.y_grid[-1] + 0.5 * dy]])
        Ix = np.digitize(np.squeeze(xc[:, 0]), x_bins) - 1
        Iy = np.digitize(np.squeeze(xc[:, 1]), y_bins) - 1
        for ii in range(Nx):
            for jj in range(Ny):
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
        if (parallel_io.rank == 0):
            range_change = np.zeros((Nx, Ny))
            for ii in range(3):
                global_displacement_sum[:, :, ii] /= global_N + 1e-10
                range_change += global_displacement_sum[:, :, ii] * self.satellite_vector[ii]

            self.displacement.append(global_displacement_sum)
            self.range_change.append(range_change)
            self.times.append(time)

    def save_csv(self, header='insar_', output_root='./results'):
        os.makedirs(output_root, exist_ok=True)

    def save_vtk(self, header='insar_', output_root='./results'):
        os.makedirs(output_root, exist_ok=True)
        print(self.times)
        print(self.displacement)
        print(self.range_change)

    def save_image(self, header='insar_', output_root='./results'):
        os.makedirs(output_root, exist_ok=True)
        fig = plot_tools.HighResPlot()
        for ii in range(len(self.times)):
            # Range change
            fig.reset()
            values = np.transpose(np.flipud(self.range_change[ii]))
            extents = [self.x_grid[0], self.x_grid[-1], self.y_grid[0], self.y_grid[-1]]
            ca = plt.imshow(values,
                            extent=extents,
                            cmap=cm.jet,
                            aspect='auto',
                            interpolation='quadric')
            plt.title('T = %1.4e (days)' % (self.times[ii] / (60 * 60 * 24)))
            plt.xlabel('X (m)')
            plt.ylabel('Y (m)')
            cb = plt.colorbar(ca)
            cb.set_label('Range Change (m)')
            fig.save('%s/%s_range_change_%03d' % (output_root, header, ii))

            # Wrapped phase
            fig.reset()
            phase = 2.0 * np.pi * values / self.satellite_wavelength
            wrapped_phase = np.angle(np.exp(1j * phase))

            extents = [self.x_grid[0], self.x_grid[-1], self.y_grid[0], self.y_grid[-1]]
            ca = plt.imshow(wrapped_phase,
                            extent=extents,
                            cmap=cm.jet,
                            aspect='auto',
                            interpolation='quadric')
            plt.title('T = %1.4e (days)' % (self.times[ii] / (60 * 60 * 24)))
            plt.xlabel('X (m)')
            plt.ylabel('Y (m)')
            cb = plt.colorbar(ca)
            cb.set_label('Phase (radians)')
            fig.save('%s/%s_wrapped_phase_%03d' % (output_root, header, ii))


