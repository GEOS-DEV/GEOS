import os
import re
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-GUI backend for cluster environments
import matplotlib.pyplot as plt
import pandas as pd
import h5py
from pathlib import Path

# Material properties
alpha_rock = 2e-6  # Thermal expansion coefficient

# Mohr-Coulomb failure criterion parameters (low strength)
cohesion_low = 3.44e6  # Pa
phi_low = 53.06  # degrees

# Mohr-Coulomb failure criterion parameters (high strength)
cohesion_high = 10.2e6  # Pa
phi_high = 39.6  # degrees

# Tolerance for finding data at specific locations
tolAngle = 2.5  # degrees
tolRadius = 0.02  # meters

# Configure matplotlib to use LaTeX rendering
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
    "mathtext.fontset": "custom",
    "mathtext.default": "regular",
})

def find_file_with_keywords(directory, *keywords):
    """Find HDF5 file in directory matching all given keywords."""
    for file in Path(directory).rglob('*'):
        if file.is_file() and all(kw in file.name for kw in keywords):
            return file
    return None

def get_data_location(data, target_radius, target_angle):
    """Extract data at specific cylindrical coordinates (radius, angle) below z=0."""
    # Filter for data below z=0
    data_targetZ = data[data["Z"] < 0.0]

    # Calculate angle in cylindrical coordinates
    data_targetZ["Angle"] = np.degrees(np.arctan2(data_targetZ['Y'], data_targetZ['X']))
    data_targetZ_sorted = data_targetZ.sort_values(by='Angle').reset_index(drop=True)

    # Filter by target angle
    data_targetAngle = data_targetZ_sorted[
        (data_targetZ_sorted["Angle"] > (target_angle - tolAngle)) &
        (data_targetZ_sorted["Angle"] < target_angle)]

    # Calculate radius
    data_targetAngle["Radius"] = np.sqrt(data_targetAngle["X"]**2 + data_targetAngle["Y"]**2)
    data_targetAngle_sorted = data_targetAngle.sort_values(by='Radius').reset_index(drop=True)

    # Filter by target radius
    data_targetRadius = data_targetAngle_sorted[
        (data_targetAngle_sorted["Radius"] > (target_radius - tolRadius)) &
        (data_targetAngle_sorted["Radius"] < (target_radius + tolRadius))]

    return data_targetRadius

def principal_components(row):
    """Calculate principal stresses and directions from 2D stress tensor."""
    stress_tensor = np.array([
        [row['compXX'], row['compXY']],
        [row['compXY'], row['compYY']]
    ])
    eigenvalues, eigenvectors = np.linalg.eigh(stress_tensor)
    return pd.Series({
        'Stress 1': eigenvalues[0],
        'Stress 2': eigenvalues[1],
        'Direction 1': eigenvectors[:, 0],
        'Direction 2': eigenvectors[:, 1]
    })


def calculate_principalStressAndDirection(data, cellCenter):
    """Calculate principal stresses from element stress data (averaged over 8 nodes)."""
    # Average stress components over 8 nodes per element
    data_xx = (data[:,0] + data[:,6] + data[:,12] + data[:,18] +
               data[:,24] + data[:,30] + data[:,36] + data[:,42]).flatten() / 8.0
    data_yy = (data[:,1] + data[:,7] + data[:,13] + data[:,19] +
               data[:,25] + data[:,31] + data[:,37] + data[:,43]).flatten() / 8.0
    data_xy = (data[:,5] + data[:,11] + data[:,17] + data[:,23] +
               data[:,29] + data[:,35] + data[:,41] + data[:,47]).flatten() / 8.0

    # Create dataframe with stress components
    data_df = pd.DataFrame({
        "compXX": data_xx,
        "compXY": data_xy,
        "compYY": data_yy
    })

    # Calculate principal stresses and directions
    data_df = data_df.apply(principal_components, axis=1)

    # Add spatial coordinates
    data_df["X"] = cellCenter[:,0].flatten()
    data_df["Y"] = cellCenter[:,1].flatten()
    data_df["Z"] = cellCenter[:,2].flatten()

    return data_df

def read_data(directory, time_option="start"):
    """Read and process stress, temperature, and bulk modulus data from HDF5 files."""
    # Read stress data
    file_rock_stress = h5py.File(find_file_with_keywords(directory, "stress", "rock"), "r")
    stress_rock = file_rock_stress["rockSolid_stress"][()]
    cellCenter_rock_data = file_rock_stress["rockSolid_stress elementCenter"][()][0]
    time = file_rock_stress["rockSolid_stress Time"][()].flatten()

    # Select timestep - find last non-zero timestep for "end" option
    if time_option == "start":
        timestep = 0
    else:
        # Find last timestep with non-zero data
        max_per_timestep = np.array([np.abs(stress_rock[i]).max() for i in range(stress_rock.shape[0])])
        non_zero_timesteps = np.where(max_per_timestep > 0)[0]
        if len(non_zero_timesteps) == 0:
            raise ValueError("No non-zero stress data found in file")
        timestep = non_zero_timesteps[-1]
        print(f"Using timestep {timestep} (last non-zero timestep) out of {len(time)} total timesteps")

    stress_rock_timestep = stress_rock[timestep, :, :]

    # Read temperature data
    file_rock_temp = h5py.File(find_file_with_keywords(directory, "temperature", "rock"), "r")
    temp_rock = file_rock_temp["temperature"][()]
    temp_rock_timestep = temp_rock[timestep, :].flatten()

    # Read bulk modulus data (use same timestep as stress)
    file_rock_bulkMod = h5py.File(find_file_with_keywords(directory, "bulkMod", "rock"), "r")
    bulkMod_rock = file_rock_bulkMod["rockSolid_bulkModulus"][()]
    bulkMod_rock_timestep = bulkMod_rock[timestep, :].flatten()

    # Calculate principal stresses and directions
    stress_rock_timestep = calculate_principalStressAndDirection(stress_rock_timestep, cellCenter_rock_data)

    # Calculate effective stresses (accounting for thermal effects)
    thermal_stress = 3 * alpha_rock * bulkMod_rock_timestep * temp_rock_timestep
    stress_rock_timestep["Effective stress 1"] = stress_rock_timestep["Stress 1"] - thermal_stress
    stress_rock_timestep["Effective stress 2"] = stress_rock_timestep["Stress 2"] - thermal_stress
    stress_rock_timestep["Temperature"] = temp_rock_timestep

    # Calculate Coulomb failure stress for low friction parameters
    stress_rock_timestep["Sigma_N low"] = 0.5 * (
        stress_rock_timestep["Effective stress 1"] + stress_rock_timestep["Effective stress 2"] +
        np.abs(stress_rock_timestep["Effective stress 1"] - stress_rock_timestep["Effective stress 2"]) *
        np.sin(np.radians(phi_low))
    )
    stress_rock_timestep["Tau low"] = 0.5 * np.abs(
        stress_rock_timestep["Effective stress 1"] - stress_rock_timestep["Effective stress 2"]
    ) * np.cos(np.radians(phi_low))
    stress_rock_timestep["Coulomb stress low"] = (
        stress_rock_timestep["Tau low"] +
        stress_rock_timestep["Sigma_N low"] * np.tan(np.radians(phi_low)) -
        cohesion_low
    )

    # Calculate Coulomb failure stress for high friction parameters
    stress_rock_timestep["Sigma_N high"] = 0.5 * (
        stress_rock_timestep["Effective stress 1"] + stress_rock_timestep["Effective stress 2"] +
        np.abs(stress_rock_timestep["Effective stress 1"] - stress_rock_timestep["Effective stress 2"]) *
        np.sin(np.radians(phi_high))
    )
    stress_rock_timestep["Tau high"] = 0.5 * np.abs(
        stress_rock_timestep["Effective stress 1"] - stress_rock_timestep["Effective stress 2"]
    ) * np.cos(np.radians(phi_high))
    stress_rock_timestep["Coulomb stress high"] = (
        stress_rock_timestep["Tau high"] +
        stress_rock_timestep["Sigma_N high"] * np.tan(np.radians(phi_high)) -
        cohesion_high
    )

    return stress_rock_timestep

def plot_mohrCircle_max(effStress_1, effStress_2, directory, cellCenter, MCCriterion):
    """Plot Mohr circle at maximum Coulomb stress location with dynamic axes."""
    fig, ax1 = plt.subplots(figsize=(8, 8))

    ax1.set_xlabel('Effective normal stress (MPa)', fontsize=12)
    ax1.set_ylabel('Shear stress (MPa)', fontsize=12)

    # Calculate Mohr circle parameters (convert to MPa)
    center = -0.5 * (effStress_1 + effStress_2) / 1e6
    radius = np.abs(effStress_1 - effStress_2) * 0.5 / 1e6

    # Plot Mohr circle
    theta = np.linspace(0, 2 * np.pi, 500)
    x = center + radius * np.cos(theta)
    y = radius * np.sin(theta)
    ax1.plot(x, y, color='r', linewidth=2, label='Mohr circle')

    # Select Mohr-Coulomb parameters
    if MCCriterion == "low":
        phi, cohesion = phi_low, cohesion_low
    elif MCCriterion == "high":
        phi, cohesion = phi_high, cohesion_high
    else:
        raise ValueError("MCCriterion must be 'low' or 'high'")

    # Determine plot range based on Mohr circle
    x_min = -effStress_2 / 1e6
    x_max = -effStress_1 / 1e6
    margin = 0.15 * (x_max - x_min)  # 15% margin
    x_plot_min = max(0, x_min - margin)
    x_plot_max = x_max + margin

    # Plot Mohr-Coulomb failure envelope
    x_mc = np.linspace(x_plot_min, x_plot_max, 100)
    y_mc = np.tan(np.radians(phi)) * x_mc + cohesion / 1e6
    ax1.plot(x_mc, y_mc, color='k', linewidth=2, label='Failure envelope')

    # Set axis limits
    y_max = max(radius * 1.2, np.tan(np.radians(phi)) * x_plot_max + cohesion / 1e6)
    ax1.set_xlim(x_plot_min, x_plot_max)
    ax1.set_ylim(0, y_max)
    ax1.set_aspect('equal', adjustable='box')
    ax1.set_title(f"Stress state at ({cellCenter[0]:.2f}, {cellCenter[1]:.2f}, {cellCenter[2]:.2f})")
    ax1.legend()
    ax1.grid(True, which='both', linestyle=':', linewidth=0.6, color='gray')

    plt.savefig(f'{directory}/mohrCircle_CoulombStress_maxLocation_{MCCriterion}.pdf',
                format='pdf', bbox_inches='tight')
    plt.close()

def plot_mohrCircle(effStress_1, effStress_2, directory, cellCenter, MCCriterion, targetRadius, targetAngle):
    """Plot Mohr circle at specific location with auto-scaled axes."""
    fig, ax1 = plt.subplots(figsize=(8, 8))

    ax1.set_xlabel('Effective normal stress (MPa)', fontsize=12)
    ax1.set_ylabel('Shear stress (MPa)', fontsize=12)

    # Calculate Mohr circle parameters (convert to MPa)
    center = -0.5 * (effStress_1 + effStress_2) / 1e6
    radius = np.abs(effStress_1 - effStress_2) * 0.5 / 1e6

    # Plot Mohr circle
    theta = np.linspace(0, 2 * np.pi, 500)
    x = center + radius * np.cos(theta)
    y = radius * np.sin(theta)
    ax1.plot(x, y, color='r', linewidth=2, label='Mohr circle')

    # Select Mohr-Coulomb parameters
    if MCCriterion == "low":
        phi, cohesion = phi_low, cohesion_low
    elif MCCriterion == "high":
        phi, cohesion = phi_high, cohesion_high
    else:
        raise ValueError("MCCriterion must be 'low' or 'high'")

    # Determine plot range based on Mohr circle
    x_min = -effStress_2 / 1e6
    x_max = -effStress_1 / 1e6
    margin = 0.15 * (x_max - x_min)  # 15% margin
    x_plot_min = max(0, x_min - margin)
    x_plot_max = x_max + margin

    # Plot Mohr-Coulomb failure envelope
    x_mc = np.linspace(x_plot_min, x_plot_max, 100)
    y_mc = np.tan(np.radians(phi)) * x_mc + cohesion / 1e6
    ax1.plot(x_mc, y_mc, color='k', linewidth=2, label='Failure envelope')

    # Set axis limits
    y_max = max(radius * 1.2, np.tan(np.radians(phi)) * x_plot_max + cohesion / 1e6)
    ax1.set_xlim(x_plot_min, x_plot_max)
    ax1.set_ylim(0, y_max)
    ax1.set_aspect('equal', adjustable='box')
    ax1.set_title(f"Stress state at ({cellCenter[0]:.2f}, {cellCenter[1]:.2f}, {cellCenter[2]:.2f})")
    ax1.legend()
    ax1.grid(True, which='both', linestyle=':', linewidth=0.6, color='gray')

    plt.savefig(f'{directory}/mohrCircle_CoulombStress_r={targetRadius}_theta={targetAngle}_{MCCriterion}.pdf',
                format='pdf', bbox_inches='tight')
    plt.close()

def write_file(effStress_1, effStress_2, T_value, targetRadius, targetAngle):
    """Append stress data to CSV summary file."""
    stress_data = pd.DataFrame({
        'EffStress 1': [effStress_1],
        'EffStress 2': [effStress_2],
        'Cooling Temp': [T_value]
    })

    file_path = f'mohrCircle_summary_varyT_r={targetRadius}_theta={targetAngle}.csv'

    if os.path.exists(file_path):
        existing_data = pd.read_csv(file_path)
        updated_data = pd.concat([existing_data, stress_data], ignore_index=True)
    else:
        updated_data = stress_data

    updated_data.to_csv(file_path, index=False) 

def main():
    """Main function to process stress data and create Mohr circle plots."""
    parser = argparse.ArgumentParser(description='Plot Mohr circles from geomechanical simulation data')
    parser.add_argument("-dir", "--target-directory", required=True,
                        help="Directory containing HDF5 output files")
    parser.add_argument("-option", "--option", required=True,
                        choices=["max", "location", "initial"],
                        help="Plot option: 'max' for maximum Coulomb stress points, "
                             "'location' for specific radius/angle, 'initial' for initial state")
    parser.add_argument("-radius", "--radius", type=float, default=None,
                        help="Target radius (required for 'location' and 'initial' options)")
    parser.add_argument("-angle", "--angle", type=float, default=None,
                        help="Target angle in degrees (required for 'location' and 'initial' options)")
    args = parser.parse_args()

    # Extract temperature from directory name
    match = re.search(r"coolT=(\d+)degC", args.target_directory)
    T_value = float(match.group(1)) if match else 0

    print("Reading data...")
    data_all = read_data(args.target_directory, "end")
    print("Stress data processed successfully")

    if args.option == "max":
        # Plot at maximum Coulomb stress locations
        print("Plotting Mohr circles at most critical points...")

        # Low friction case
        idx_low = data_all['Coulomb stress low'].idxmax()
        cellCenter_max_low = np.array([
            data_all.loc[idx_low, 'X'],
            data_all.loc[idx_low, 'Y'],
            data_all.loc[idx_low, 'Z']
        ])
        effStress_1_low = data_all.loc[idx_low, 'Effective stress 1']
        effStress_2_low = data_all.loc[idx_low, 'Effective stress 2']
        plot_mohrCircle_max(effStress_1_low, effStress_2_low, args.target_directory,
                           cellCenter_max_low, "low")

        # High friction case
        idx_high = data_all['Coulomb stress high'].idxmax()
        cellCenter_max_high = np.array([
            data_all.loc[idx_high, 'X'],
            data_all.loc[idx_high, 'Y'],
            data_all.loc[idx_high, 'Z']
        ])
        effStress_1_high = data_all.loc[idx_high, 'Effective stress 1']
        effStress_2_high = data_all.loc[idx_high, 'Effective stress 2']
        plot_mohrCircle_max(effStress_1_high, effStress_2_high, args.target_directory,
                           cellCenter_max_high, "high")

    elif args.option == "location":
        # Plot at specific location
        if args.radius is None or args.angle is None:
            raise ValueError("Radius and angle are required for 'location' option")

        print(f"Getting stress at r={args.radius}, theta={args.angle}...")
        data_loc = get_data_location(data_all, args.radius, args.angle)

        cellCenter = np.array([data_loc["X"].iloc[0], data_loc["Y"].iloc[0], data_loc["Z"].iloc[0]])
        effStress_1 = data_loc['Effective stress 1'].iloc[0]
        effStress_2 = data_loc['Effective stress 2'].iloc[0]

        print("Plotting Mohr circles...")
        plot_mohrCircle(effStress_1, effStress_2, args.target_directory, cellCenter,
                       "low", args.radius, args.angle)
        plot_mohrCircle(effStress_1, effStress_2, args.target_directory, cellCenter,
                       "high", args.radius, args.angle)

        write_file(effStress_1, effStress_2, T_value, args.radius, args.angle)

    elif args.option == "initial":
        # Get initial stress state
        if args.radius is None or args.angle is None:
            raise ValueError("Radius and angle are required for 'initial' option")

        print("Reading initial stress state...")
        data_all_initial = read_data(args.target_directory, "start")

        print(f"Getting initial stress at r={args.radius}, theta={args.angle}...")
        data_loc = get_data_location(data_all_initial, args.radius, args.angle)

        effStress_1 = data_loc['Effective stress 1'].iloc[0]
        effStress_2 = data_loc['Effective stress 2'].iloc[0]

        write_file(effStress_1, effStress_2, "Initial State", args.radius, args.angle)

    print("Done!")


if __name__ == '__main__':
    main()