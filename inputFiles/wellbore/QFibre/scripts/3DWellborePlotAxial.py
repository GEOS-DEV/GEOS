import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import h5py
from pathlib import Path

targetAngle = 0
tolAngle = 2.5
bulkMod_casing = 166.66666667e9
bulkMod_cement = 6.88888889e9
bulkMod_seal = 13.33333333e9
alpha_casing = 1.2e-5
alpha_cement = 2e-5
alpha_seal = 2e-5

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

def get_data_axial(cellCenter, data, targetRadius, tolRadius):
    data_df = pd.DataFrame({})
    data_df["Variable"] = data
    data_df["X"] = cellCenter[:,0].flatten() 
    data_df["Y"] = cellCenter[:,1].flatten() 
    data_df["Z"] = cellCenter[:,2].flatten() 
    data_df["Radius"] = np.sqrt(data_df["X"]**2 + data_df["Y"]**2)
    data_df_targetR = data_df[
        (data_df["Radius"] > targetRadius) &
        (data_df["Radius"] < (targetRadius + tolRadius))
    ]
    data_df_targetR["Angle"] = np.degrees(np.arctan2(data_df_targetR['Y'], data_df_targetR['X']))
    data_df_targetR_sorted = data_df_targetR.sort_values(by='Angle').reset_index(drop=True)

    data_df_targetR_targetAngle = data_df_targetR_sorted[
        (data_df_targetR_sorted["Angle"] > (targetAngle - tolAngle)) &
        (data_df_targetR_sorted["Angle"] < targetAngle)].reset_index(drop=True)
    
    data_df_targetR_targetAngle_sorted = data_df_targetR_targetAngle.sort_values(by='Z').reset_index(drop=True)

    return data_df_targetR_targetAngle_sorted

def read_and_plot_scalar( directory, variable ):
    file_casing = h5py.File(find_file_with_keywords(directory, variable, "casing"), "r")
    file_cement = h5py.File(find_file_with_keywords(directory, variable, "cement"), "r")
    file_seal = h5py.File(find_file_with_keywords(directory, variable, "seal"), "r")

    # Process for casing
    #-------------------
    data_casing = file_casing[variable][()]
    cellCenter_casing_data = file_casing[variable+" elementCenter"][()][0]
    time_data = file_casing[variable+" Time"][()].flatten()

    # We pick two time instances
    indice1, indice2 = int(len(time_data)/2), len(time_data)
    time1, time2 = time_data[indice1-1], time_data[indice2-1]

    data_casing_time1, data_casing_time2 = data_casing[indice1-1,:].flatten(), data_casing[indice2-1,:].flatten()

    axialData_casing_time1 = get_data_axial(cellCenter_casing_data, data_casing_time1, 0.1, 0.003)
    axialData_casing_time2 = get_data_axial(cellCenter_casing_data, data_casing_time2, 0.1, 0.003)

    # Process for cement
    #-------------------
    data_cement = file_cement[variable][()]
    cellCenter_cement_data = file_cement[variable+" elementCenter"][()][0]

    data_cement_time1, data_cement_time2 = data_cement[indice1-1,-25600:].flatten(), data_cement[indice2-1,-25600:].flatten()

    axialData_cement_time1 = get_data_axial(cellCenter_cement_data[-25600:], data_cement_time1, 0.12, 0.006)
    axialData_cement_time2 = get_data_axial(cellCenter_cement_data[-25600:], data_cement_time2, 0.12, 0.006)

    # Process for seal
    #-------------------
    data_seal = file_seal[variable][()]
    cellCenter_seal_data = file_seal[variable+" elementCenter"][()][0]

    data_seal_time1, data_seal_time2 = data_seal[indice1-1,:].flatten(), data_seal[indice2-1,:].flatten()

    axialData_seal_time1 = get_data_axial(cellCenter_seal_data, data_seal_time1, 0.16, 0.015)
    axialData_seal_time2 = get_data_axial(cellCenter_seal_data, data_seal_time2, 0.16, 0.015)

    plot_all(axialData_casing_time1, axialData_cement_time1, axialData_seal_time1, time1, variable, directory)
    plot_all(axialData_casing_time2, axialData_cement_time2, axialData_seal_time2, time2, variable, directory)

    return True 

def read_and_plot_tensor( directory, variable ):
    # For averageStrain, the files are named with "strain" prefix
    file_keyword = "strain" if variable == "averageStrain" else variable
    file_casing = h5py.File(find_file_with_keywords(directory, file_keyword, "casing"), "r")
    file_cement = h5py.File(find_file_with_keywords(directory, file_keyword, "cement"), "r")
    file_seal = h5py.File(find_file_with_keywords(directory, file_keyword, "seal"), "r")

    # Process for casing
    #-------------------
    variableName_casing = variable 
    if variable == "stress":
        variableName_casing = "casingSolid_"+variable
    data_casing = file_casing[variableName_casing][()]
    cellCenter_casing_data = file_casing[variableName_casing+" elementCenter"][()][0]
    time_data = file_casing[variableName_casing+" Time"][()].flatten()

    # We pick two time instances
    indice1, indice2 = int(len(time_data)/2), len(time_data)
    time1, time2 = time_data[indice1-1], time_data[indice2-1]

    data_casing_time1, data_casing_time2 = data_casing[indice1-1,:,2].flatten(), data_casing[indice2-1,:,2].flatten()

    axialData_casing_time1 = get_data_axial(cellCenter_casing_data, data_casing_time1, 0.1, 0.003)
    axialData_casing_time2 = get_data_axial(cellCenter_casing_data, data_casing_time2, 0.1, 0.003)

    if variable == "stress":
        file_casing_temp = h5py.File(find_file_with_keywords(directory, "temperature", "casing"), "r")
        temp_casing = file_casing_temp["temperature"][()]
        temp_casing_time1, temp_casing_time2 = temp_casing[indice1-1,:].flatten(), temp_casing[indice2-1,:].flatten()
        axialTemp_casing_time1 = get_data_axial(cellCenter_casing_data, temp_casing_time1, 0.1, 0.003)
        axialTemp_casing_time2 = get_data_axial(cellCenter_casing_data, temp_casing_time2, 0.1, 0.003)

        axialData_casing_time1["EffAxialStress"] = axialData_casing_time1["Variable"] - 3*alpha_casing*bulkMod_casing*axialTemp_casing_time1["Variable"]
        axialData_casing_time2["EffAxialStress"] = axialData_casing_time2["Variable"] - 3*alpha_casing*bulkMod_casing*axialTemp_casing_time2["Variable"]

    # Process for cement
    #-------------------
    variableName_cement = variable 
    if variable == "stress":
        variableName_cement = "cementSolid_"+variable
    data_cement = file_cement[variableName_cement][()]
    cellCenter_cement_data = file_cement[variableName_cement+" elementCenter"][()][0]

    data_cement_time1, data_cement_time2 = data_cement[indice1-1,-25600:,2].flatten(), data_cement[indice2-1,-25600:,2].flatten()

    axialData_cement_time1 = get_data_axial(cellCenter_cement_data[-25600:], data_cement_time1, 0.12, 0.006)
    axialData_cement_time2 = get_data_axial(cellCenter_cement_data[-25600:], data_cement_time2, 0.12, 0.006)

    if variable == "stress":
        file_cement_temp = h5py.File(find_file_with_keywords(directory, "temperature", "cement"), "r")
        temp_cement = file_cement_temp["temperature"][()]
        temp_cement_time1, temp_cement_time2 = temp_cement[indice1-1,-25600:].flatten(), temp_cement[indice2-1,-25600:].flatten()
        axialTemp_cement_time1 = get_data_axial(cellCenter_cement_data[-25600:], temp_cement_time1, 0.12, 0.006)
        axialTemp_cement_time2 = get_data_axial(cellCenter_cement_data[-25600:], temp_cement_time2, 0.12, 0.006)

        axialData_cement_time1["EffAxialStress"] = axialData_cement_time1["Variable"] - 3*alpha_cement*bulkMod_cement*axialTemp_cement_time1["Variable"]
        axialData_cement_time2["EffAxialStress"] = axialData_cement_time2["Variable"] - 3*alpha_cement*bulkMod_cement*axialTemp_cement_time2["Variable"]

    # Process for seal
    #-------------------
    variableName_seal = variable 
    if variable == "stress":
        variableName_seal = "sealSolid_"+variable
    data_seal = file_seal[variableName_seal][()]
    cellCenter_seal_data = file_seal[variableName_seal+" elementCenter"][()][0]

    data_seal_time1, data_seal_time2 = data_seal[indice1-1,:,2].flatten(), data_seal[indice2-1,:,2].flatten()

    axialData_seal_time1 = get_data_axial(cellCenter_seal_data, data_seal_time1, 0.16, 0.015)
    axialData_seal_time2 = get_data_axial(cellCenter_seal_data, data_seal_time2, 0.16, 0.015)

    if variable == "stress":
        file_seal_temp = h5py.File(find_file_with_keywords(directory, "temperature", "seal"), "r")
        temp_seal = file_seal_temp["temperature"][()]
        temp_seal_time1, temp_seal_time2 = temp_seal[indice1-1,:].flatten(), temp_seal[indice2-1,:].flatten()
        axialTemp_seal_time1 = get_data_axial(cellCenter_seal_data, temp_seal_time1, 0.16, 0.015)
        axialTemp_seal_time2 = get_data_axial(cellCenter_seal_data, temp_seal_time2, 0.16, 0.015)

        axialData_seal_time1["EffAxialStress"] = axialData_seal_time1["Variable"] - 3*alpha_seal*bulkMod_seal*axialTemp_seal_time1["Variable"]
        axialData_seal_time2["EffAxialStress"] = axialData_seal_time2["Variable"] - 3*alpha_seal*bulkMod_seal*axialTemp_seal_time2["Variable"]
    
    plot_all(axialData_casing_time1, axialData_cement_time1, axialData_seal_time1, time1, variable, directory)
    plot_all(axialData_casing_time2, axialData_cement_time2, axialData_seal_time2, time2, variable, directory)

    return True

def plot_all(casing_data, cement_data, seal_data, time, variable, directory):
    fig, ax1 = plt.subplots()
    
    ax1.set_ylabel('Depth (m)', fontsize=12)
    if variable=="temperature":
        ax1.plot(casing_data['Variable']-273, casing_data['Z'], 's', color='r', markerfacecolor='none', label='Inside casing')
        ax1.plot(cement_data['Variable']-273, cement_data['Z'], 's', color='b', markerfacecolor='none', label='Behind casing')
        ax1.plot(seal_data['Variable']-273, seal_data['Z'], 's', color='g', markerfacecolor='none', label='Inside formation')
        ax1.set_xlabel('Temperature ($^\\circ$C)', fontsize=12)
    elif variable=="pressure":
        # ax1.plot(casing_data['Variable']/1e6, casing_data['Z'], 's', color='r', markerfacecolor='none', label='Inside casing')
        ax1.plot(cement_data['Variable']/1e6, cement_data['Z'], 's', color='b', markerfacecolor='none', label='Behind casing')
        ax1.plot(seal_data['Variable']/1e6, seal_data['Z'], 's', color='g', markerfacecolor='none', label='Inside formation')
        ax1.set_xlabel('Pressure (MPa)', fontsize=12)
    elif variable=="averageStrain":
        ax1.plot(casing_data['Variable']*100, casing_data['Z'],  's', color='r', markerfacecolor='none', label='Inside casing')
        ax1.plot(cement_data['Variable']*100, cement_data['Z'], 's', color='b', markerfacecolor='none', label='Behind casing')
        ax1.plot(seal_data['Variable']*100, seal_data['Z'], 's', color='g', markerfacecolor='none', label='Inside formation')
        ax1.set_xlabel('Axial strain (\%)', fontsize=12)
    elif variable=="stress":
        ax1.plot(casing_data['EffAxialStress']/1e6, casing_data['Z'],  's', color='r', markerfacecolor='none', label='Inside casing')
        ax1.plot(cement_data['EffAxialStress']/1e6, cement_data['Z'], 's', color='b', markerfacecolor='none', label='Behind casing')
        ax1.plot(seal_data['EffAxialStress']/1e6, seal_data['Z'], 's', color='g', markerfacecolor='none', label='Inside formation')
        ax1.set_xlabel('Effective axial stress (MPa)', fontsize=12)

    ax1.legend(frameon=False, fontsize=12, loc='best')
    fig.tight_layout()
    plt.savefig(directory+'/'+variable+'_axialDistribution_'+str(time)+'s.pdf', format='pdf', bbox_inches='tight')

    return True

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-dir", "--target-directory", help=" target directory ")
    parser.add_argument("-var", "--target-variable", help=" target variable ")
    args = parser.parse_args()

    if args.target_variable == "stress" or args.target_variable == "averageStrain":
        read_and_plot_tensor( args.target_directory, args.target_variable )
    elif args.target_variable == "temperature" or args.target_variable == "pressure":
        read_and_plot_scalar( args.target_directory, args.target_variable )
    else:
        assert False, "No such variable supported"
    
if __name__ == '__main__':
    main()