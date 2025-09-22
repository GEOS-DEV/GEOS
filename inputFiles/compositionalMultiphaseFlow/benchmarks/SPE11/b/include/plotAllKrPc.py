import matplotlib.pyplot as plt
import numpy as np
import argparse

def main():
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Script to generate figure from tutorial.")

    # Add arguments to accept individual file paths
    parser.add_argument('--outputDir', help='Path to output directory', default='..')

    # Parse the command-line arguments
    args = parser.parse_args()

    # File path
    base = args.outputDir

    fsize = 30
    msize = 12
    lw = 4
    fig, ax = plt.subplots(1, 2, figsize=(32, 18))
    cmap = plt.get_cmap("tab10")

    for i in range(1, 7):
        s = np.loadtxt(f'{base}/tables/gasSaturation_spe11b_{i}.txt')
        sn = np.loadtxt(f'{base}/tables/waterSaturation_spe11b_{i}.txt')
        spc = np.loadtxt(f'{base}/tables/waterPCSaturation_spe11b_{i}.txt')
        kr = np.loadtxt(f'{base}/tables/gasRelperm_spe11b_{i}.txt')
        krn = np.loadtxt(f'{base}/tables/waterRelperm_spe11b_{i}.txt')
        pc = np.loadtxt(f'{base}/tables/waterCapPres_spe11b_{i}.txt')
        if not i == 5:
            ax[0].plot(s, kr, color=cmap(-1), label=f'facies {i}', lw=lw)
            ax[0].plot(sn[::-1], krn, color=cmap(-1), lw=lw, ls=':')
            ax[1].plot(spc, pc, color=cmap(-1), label=f'facies {i}', lw=lw)
        else:
            ax[0].plot(s, kr, label=f'facies {i}', lw=lw, ls='-', color='red')
            ax[0].plot(sn[::-1], krn, lw=lw, ls='-',color='red')
            ax[1].plot(spc, pc, label=f'facies {i}', lw=lw, ls='-',color='red')

    ax[0].set_xlabel('Relative Permeability [-]', size=fsize, weight="bold")
    ax[0].set_ylabel('Wetting saturation [-]', size=fsize, weight="bold")
    ax[0].legend(loc='upper left')
    ax[1].set_xlabel('Capillary pressure [Pa]', size=fsize, weight="bold")
    ax[1].set_ylabel('Wetting saturation [-]', size=fsize, weight="bold")
    ax[1].legend(loc='upper left')
    plt.show()


if __name__ == "__main__":
    main()
