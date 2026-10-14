import h5py
import numpy
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as patches
from matplotlib.ticker import MaxNLocator
import argparse

FONT_SIZE = 22

def read_dataset_and_grid(hdf_file: str, name: str):
    """Reads dataset and grid properties assuming attributes always exist."""
    with h5py.File(hdf_file, 'r') as h5f:
        if name not in h5f:
            raise ValueError(f"Dataset '{name}' not found in {hdf_file}.")
        data = h5f[name][:]
        dx, dz = h5f.attrs['grid_spacing']
        xmin, zmin = h5f.attrs['grid_origin']
    return data, dx, dz, xmin, zmin

def setup_plot_layout(nz: int, nx: int, dx: float, dz: float, xmin: float, zmin: float):
    """Sets up the standardized figure, image axis, and layout dimensions."""
    data_aspect = (nz * dz) / (nx * dx)
    
    dpi = 100
    pad_in = 3.0 / dpi
    
    img_width_in = 12.0
    img_height_in = img_width_in * data_aspect
    
    fig_width = img_width_in + 2.0 * pad_in
    legend_height = 1.2 
    fig_height = img_height_in + legend_height + pad_in
    
    fig = plt.gcf()
    fig.clf()
    fig.set_size_inches(fig_width, fig_height)
    fig.set_facecolor('white')
    fig.set_dpi(dpi)

    ax_img_x = pad_in / fig_width
    ax_img_y = legend_height / fig_height
    ax_img_w = img_width_in / fig_width
    ax_img_h = img_height_in / fig_height
    
    ax_img = fig.add_axes([ax_img_x, ax_img_y, ax_img_w, ax_img_h])
    extent = [xmin, xmin + nx * dx, zmin, zmin + nz * dz]
    
    return fig, ax_img, extent, fig_height

def plot_attribute(hdf_file: str, component_name: str = "", units: str = "-"):
    """
    Reads the 'attribute' array and generates a discrete colored plot.
    """
    attribute, dx, dz, xmin, zmin = read_dataset_and_grid(hdf_file, 'attribute')
        
    valid_data = attribute[~numpy.isnan(attribute)]
    if len(valid_data) == 0:
        raise ValueError("The 'attribute' array contains only NaNs.")
        
    unique_vals = numpy.sort(numpy.unique(valid_data))
    target_colors = ['#00BFFF', '#FF0000', '#00FF00']
    labels = ['SAND', 'CLAY', 'CAPROCK']
    
    n_vals = len(unique_vals)
    if n_vals > len(target_colors):
        target_colors += ['#808080'] * (n_vals - len(target_colors))
        labels += [f'Type {i}' for i in range(len(target_colors), n_vals)]
        
    cmap = mcolors.ListedColormap(target_colors[:n_vals])
    cmap.set_bad(color='none')
    
    if n_vals == 1:
        boundaries = [unique_vals[0] - 0.5, unique_vals[0] + 0.5]
    else:
        diffs = numpy.diff(unique_vals)
        min_diff = numpy.min(diffs) / 2
        boundaries = [v - min_diff for v in unique_vals] + [unique_vals[-1] + min_diff]
        
    norm = mcolors.BoundaryNorm(boundaries, cmap.N)

    nz, nx = attribute.shape
    fig, ax_img, extent, fig_height = setup_plot_layout(nz, nx, dx, dz, xmin, zmin)
    
    ax_img.imshow(attribute, origin='lower', extent=extent, cmap=cmap, norm=norm, aspect='auto')
    ax_img.axis('off')
    
    ax_leg_y = 0.15 / fig_height
    ax_leg_h = 0.75 / fig_height
    ax_leg = fig.add_axes([0.3, ax_leg_y, 0.4, ax_leg_h])
    ax_leg.axis('off')
    
    title = 'Rock Type'
    if component_name:
        title += f" {component_name}"
    title += f" [{units}]"
    
    ax_leg.text(0.5, 0.9, title, ha='center', va='center', fontsize=FONT_SIZE)
    
    gap = 0.02
    box_width = (1.0 - gap * (n_vals - 1)) / n_vals
    for i, (val, color, label) in enumerate(zip(unique_vals, target_colors[:n_vals], labels[:n_vals])):
        x_pos = i * (box_width + gap)
        rect = patches.Rectangle((x_pos, 0.25), box_width, 0.4, facecolor=color, edgecolor='none')
        ax_leg.add_patch(rect)
        ax_leg.text(x_pos + box_width / 2, 0.1, label, ha='center', va='top', fontsize=FONT_SIZE)
        
    ax_leg.set_xlim(0, 1)
    ax_leg.set_ylim(0, 1)

    return fig

def plot_property(hdf_file: str, name: str, component_idx: int = 0, component_name: str = "", units: str = "-", value_range = None):
    """
    Reads a real-valued property array and generates a continuous colored plot.
    """
    data, dx, dz, xmin, zmin = read_dataset_and_grid(hdf_file, name)
    
    if data.ndim == 3:
        if component_idx >= data.shape[2]:
            raise ValueError(f"Component {component_idx} out of bounds for dataset {name}.")
        data = data[:, :, component_idx]
        
    masked_data = numpy.ma.masked_invalid(data)
    
    vmin = value_range[0] if value_range else None
    vmax = value_range[-1] if value_range else None
    
    cmap = plt.cm.coolwarm
    cmap.set_bad(color='none')
    
    nz, nx = data.shape
    fig, ax_img, extent, fig_height = setup_plot_layout(nz, nx, dx, dz, xmin, zmin)
    
    im = ax_img.imshow(masked_data, origin='lower', extent=extent, cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
    ax_img.axis('off')

    ax_cbar_y = 0.15 / fig_height
    ax_cbar_h = 0.2 / fig_height
    ax_cbar = fig.add_axes([0.3, ax_cbar_y, 0.4, ax_cbar_h])
    
    cbar = fig.colorbar(im, cax=ax_cbar, orientation='horizontal')
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.tick_params(labelsize=FONT_SIZE, length=2, width=1, direction='out')
    cbar.locator = MaxNLocator(nbins=5)
    cbar.update_ticks()
    
    title = name
    if component_name:
        title += f" {component_name}"
    title += f" [{units}]"
    
    title_y = ax_cbar_y + ax_cbar_h + (0.55 / fig_height)
    fig.text(0.5, title_y, title, ha='center', va='center', fontsize=FONT_SIZE)

    return fig

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot HDF5 datasets to match Paraview style.")
    parser.add_argument("-i", "--input", type=str, required=True, help="Path to the input HDF5 file.")
    parser.add_argument("-p", "--property", type=str, required=True, help="Name of the property to plot (e.g., 'attribute' or 'globalCompFraction').")
    parser.add_argument("-u", "--units", type=str, default="-", help="Units for the property (default: '-').")
    parser.add_argument("-n", "--number", type=int, default=0, help="Component number for multi-component datasets (default: 0).")
    parser.add_argument("-c", "--component", type=str, default="", help="Name of the component (default: '').")
    
    args = parser.parse_args()
    
    if args.property == "attribute":
        plot_attribute(args.input, component_name=args.component, units=args.units)
    else:
        plot_property(args.input, args.property, component_idx=args.number, component_name=args.component, units=args.units)
    
    plt.show()
