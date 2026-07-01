import os
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import h5py
from pathlib import Path

tolAngle = 2.5
tolRadius = 0.01

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",  # Use a serif font family (required for LaTeX rendering)
    "font.serif": ["Palatino"],  # Specify the custom LaTeX font
    "mathtext.fontset": "custom",  # Use the custom LaTeX font for math text
    "mathtext.default": "regular",  # Set the math font to regular
})

def find_file_with_keywords(directory, *keywords):
    for file in Path(directory).rglob('*'):
        if file.is_file() and all(kw in file.name for kw in keywords):
            return file
    return None

def get_data_location(data, target_radius, target_angle):
    data_df_targetZ = data[data["Z"] < 0.0].copy()
    data_df_targetZ["Angle"] = np.degrees(np.arctan2(data_df_targetZ['Y'], data_df_targetZ['X']))
    data_df_targetZ_sorted = data_df_targetZ.sort_values(by='Angle').reset_index(drop=True)

    data_df_targetZ_targetAngle = data_df_targetZ_sorted[
        (data_df_targetZ_sorted["Angle"] > (target_angle - tolAngle)) &
        (data_df_targetZ_sorted["Angle"] < target_angle)].copy()
    data_df_targetZ_targetAngle["Radius"] = np.sqrt(data_df_targetZ_targetAngle["X"]**2 + data_df_targetZ_targetAngle["Y"]**2)
    data_df_targetZ_targetAngle_sorted = data_df_targetZ_targetAngle.sort_values(by='Radius').reset_index(drop=True)

    data_df_targetZ_targetAngle_targetRadius = data_df_targetZ_targetAngle_sorted[
        (data_df_targetZ_targetAngle_sorted["Radius"] > (target_radius - tolRadius)) &
        (data_df_targetZ_targetAngle_sorted["Radius"] < (target_radius + tolRadius))]

    return data_df_targetZ_targetAngle_targetRadius

def principal_components(row):
    A = np.array([
        [row['compXX'], row['compXY']],
        [row['compXY'], row['compYY']]
    ])
    eigenvalues, eigenvectors = np.linalg.eigh(A)
    return pd.Series({
        'Effective stress 1': eigenvalues[0],
        'Effective stress 2': eigenvalues[1],
        'Direction 1': eigenvectors[:, 0],
        'Direction 2': eigenvectors[:, 1]
    })


def calculate_principalStressAndDirection(data, cellCenter):
    data_xx = (data[:,0].flatten() + data[:,6].flatten() + data[:,12].flatten() + data[:,18].flatten() + data[:,24].flatten() + data[:,30].flatten() + data[:,36].flatten() + data[:,42].flatten())/8.
    data_yy = (data[:,1].flatten() + data[:,7].flatten() + data[:,13].flatten() + data[:,19].flatten() + data[:,25].flatten() + data[:,31].flatten() + data[:,37].flatten() + data[:,43].flatten())/8.
    data_xy = (data[:,5].flatten() + data[:,11].flatten() + data[:,17].flatten() + data[:,23].flatten() + data[:,29].flatten() + data[:,35].flatten() + data[:,41].flatten() + data[:,47].flatten())/8.

    data_df = pd.DataFrame({})
    data_df["compXX"] = data_xx
    data_df["compXY"] = data_xy
    data_df["compYY"] = data_yy

    data_df = data_df.apply(principal_components, axis=1)

    data_df["X"] = cellCenter[:,0].flatten() 
    data_df["Y"] = cellCenter[:,1].flatten() 
    data_df["Z"] = cellCenter[:,2].flatten() 

    return data_df

def read_data( stress_rock, cellCenter_rock_data, timestep ):   
    stress_rock_timestep = stress_rock[timestep,:,:]

    # Get principal stresses and directions
    stress_rock_timestep = calculate_principalStressAndDirection(stress_rock_timestep, cellCenter_rock_data)

    return stress_rock_timestep

def plot_mohrCircle(effStress_1, effStress_2, timeToPlot, directory, cellCenter, targetRadius, targetAngle):
    _, ax1 = plt.subplots()

    ax1.set_xlabel('Effective normal stress (MPa)', fontsize=12)
    ax1.set_ylabel('Shear stress (MPa)', fontsize=12)

    colors = ['#08519C', '#2171B5', '#4292C6', '#6BAED6']

    for i in range(len(effStress_1)):
        center = -0.5*(effStress_1[i] + effStress_2[i])/1e6
        radius_mc = np.abs(effStress_1[i] - effStress_2[i])*0.5/1e6

        # Create circle
        theta = np.linspace(0, 2 * np.pi, 500)
        x = center + radius_mc * np.cos(theta)
        y = radius_mc * np.sin(theta)
        ax1.plot(x, y, color=colors[i], label=f"Time = {timeToPlot[i]}s")

    ax1.set_ylim(0, 20)
    ax1.set_aspect('equal', adjustable='box')
    ax1.set_title(f"Stress state at ({cellCenter[0]:.2f}, {cellCenter[1]:.2f}, {cellCenter[2]:.2f})")
    ax1.grid(True, which='both', linestyle=':', linewidth=0.6, color='gray')
    ax1.legend(frameon=True, fontsize=10, loc='best')

    plt.savefig(f'{directory}/mohrCircle_CoulombStress_r={targetRadius}_theta={targetAngle}.pdf',
                format='pdf', bbox_inches='tight')
    plt.close()

def write_file(effStress_1, effStress_2, timeToPlot, targetRadius, targetAngle):
    stress_data = pd.DataFrame({
        'EffStress 1': effStress_1,
        'EffStress 2': effStress_2,
        'Time': timeToPlot
    })

    file_path = f'mohrCircle_summary_DP_r={targetRadius}_theta={targetAngle}.csv'

    if os.path.exists(file_path):
        existing_data = pd.read_csv(file_path)
        updated_data = pd.concat([existing_data, stress_data], ignore_index=True)
    else:
        updated_data = stress_data

    updated_data.to_csv(file_path, index=False) 

def main():
    parser = argparse.ArgumentParser(description='Plot Mohr Circle from stress data')
    parser.add_argument("-dir", "--target-directory", help="Target directory containing results",
                        default="results-plastic")
    parser.add_argument("-radius", "--radius", type=float, required=True,
                        help='The value of specific radius')
    parser.add_argument("-angle", "--angle", type=float, required=True,
                        help='The value of specific angle')
    args = parser.parse_args()

    print("START READING DATA")

    # Read stress data
    stress_file = find_file_with_keywords(args.target_directory, "stress", "rock")
    if stress_file is None:
        raise FileNotFoundError(f"No stress file found in {args.target_directory}")

    with h5py.File(stress_file, "r") as file_rock_stress:
        stress_rock = file_rock_stress["rockSolid_stress"][()]
        cellCenter_rock_data = file_rock_stress["rockSolid_stress elementCenter"][()][0]
        time = file_rock_stress["rockSolid_stress Time"][()].flatten()

    data_t0 = read_data(stress_rock, cellCenter_rock_data, 0)
    data_t1 = read_data(stress_rock, cellCenter_rock_data, int(len(time)/3))
    data_t2 = read_data(stress_rock, cellCenter_rock_data, int(len(time)*2/3))
    data_tend = read_data(stress_rock, cellCenter_rock_data, len(time)-1)

    timeToPlot = np.array([time[0], time[int(len(time)/3)],
                           time[int(len(time)*2/3)], time[-1]])

    print("STRESS DATA READ AND PROCESSED")

    print("GETTING STRESS AT THE SPECIFIED LOCATION")
    data_loc_t0 = get_data_location(data_t0, args.radius, args.angle)
    data_loc_t1 = get_data_location(data_t1, args.radius, args.angle)
    data_loc_t2 = get_data_location(data_t2, args.radius, args.angle)
    data_loc_tend = get_data_location(data_tend, args.radius, args.angle)

    print(data_loc_t0)
    print(data_loc_tend)

    cellCenter = np.array([data_loc_t0["X"].iloc[0], data_loc_t0["Y"].iloc[0],
                          data_loc_t0["Z"].iloc[0]])
    effStress_1 = np.array([
        data_loc_t0['Effective stress 1'].iloc[0],
        data_loc_t1['Effective stress 1'].iloc[0],
        data_loc_t2['Effective stress 1'].iloc[0],
        data_loc_tend['Effective stress 1'].iloc[0]
    ])

    effStress_2 = np.array([
        data_loc_t0['Effective stress 2'].iloc[0],
        data_loc_t1['Effective stress 2'].iloc[0],
        data_loc_t2['Effective stress 2'].iloc[0],
        data_loc_tend['Effective stress 2'].iloc[0]
    ])

    print("PLOTTING MOHR CIRCLE AT THE SPECIFIED LOCATION")
    plot_mohrCircle(effStress_1, effStress_2, timeToPlot, args.target_directory,
                   cellCenter, args.radius, args.angle)

    write_file(effStress_1, effStress_2, timeToPlot, args.radius, args.angle)

    print("DONE")


if __name__ == '__main__':
    main()