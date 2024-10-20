#!/usr/bin/python

import matplotlib.pyplot as plt
import seaborn as sns

# Set the aesthetic style of the plots
sns.set(style="white")

# Raw data extracted from the following runs:
#   Mechanics:  frontier-[2002292-2002296]
#   Comp. flow: frontier-[2002210-2002213]
def weak_scaling_plot(savename, lw=3, mks=11):
    fig, axs = plt.subplots(1, 2, figsize=[18, 6])  # Create a figure and two subplots (rows)

    # Plot 1 (Mechanics)
    ranks1 = [8, 64, 512, 4096, 32768]
    global_dofs1 = ['20.7M', '162M', '1.3B', '10.2B', '81.3B']
    geos_times1 = [1.7805, 1.9492, 1.8850, 1.9575, 2.3137]
    matrix_times1 = [0.4270, 0.4742, 0.4755, 0.4963, 0.7132]
    hypre_setup_times1 = [0.9412, 1.1372, 1.3507, 1.6657, 2.2864]
    hypre_solve_times1 = [0.6751, 1.0381, 1.4013, 2.0558, 3.0379]
    total_times1 = [a + b + c + d for a, b, c, d in zip(geos_times1, matrix_times1, hypre_setup_times1, hypre_solve_times1)]

    axs[0].plot(ranks1, geos_times1, 'o-', label='GEOS', linewidth=lw, markersize=mks)
    axs[0].plot(ranks1, matrix_times1, 's-', label='Matrix creation', linewidth=lw, markersize=mks)
    axs[0].plot(ranks1, hypre_setup_times1, '^-', label='Hypre setup', linewidth=lw, markersize=mks)
    axs[0].plot(ranks1, hypre_solve_times1, 'v-', label='Hypre solve', linewidth=lw, markersize=mks)
    axs[0].plot(ranks1, total_times1, 'd-', label='Total', linewidth=lw, markersize=mks)
    axs[0].set_title('(a) Mechanics', fontsize=23, pad=10)
    axs[0].set_xlabel('Number of GPUs (Global DOFs)', fontsize=20)
    axs[0].set_ylabel('Time [s]', fontsize=19)
    axs[0].set_xscale('log', base=2)
    axs[0].set_xticks(ranks1)
    axs[0].set_xticklabels([f"{r:,} ({d})" for r, d in zip(ranks1, global_dofs1)], fontsize=17, rotation=30)
    axs[0].tick_params(axis='y', labelsize=17)
    axs[0].set_ylim(bottom=0)
    axs[0].grid(True, which='both', linestyle='--')

    # Plot 2 (Compositional flow)
    ranks2 = [4, 32, 256, 2048]
    global_dofs2 = ['19.2M', '153M', '1.2B', '9.8B']
    geos_times2 = [0.1909, 0.2051, 0.2109, 0.2180]
    matrix_times2 = [0.1544, 0.1638, 0.1672, 0.1693]
    hypre_setup_times2 = [0.5262, 0.6353, 0.7320, 0.8254]
    hypre_solve_times2 = [0.1053, 0.1253, 0.1487, 0.1797]
    total_times2 = [a + b + c + d for a, b, c, d in zip(geos_times2, matrix_times2, hypre_setup_times2, hypre_solve_times2)]

    axs[1].plot(ranks2, geos_times2, 'o-', label='GEOS', linewidth=lw, markersize=mks)
    axs[1].plot(ranks2, matrix_times2, 's-', label='Matrix creation', linewidth=lw, markersize=mks)
    axs[1].plot(ranks2, hypre_setup_times2, '^-', label='Hypre setup', linewidth=lw, markersize=mks)
    axs[1].plot(ranks2, hypre_solve_times2, 'v-', label='Hypre solve', linewidth=lw, markersize=mks)
    axs[1].plot(ranks2, total_times2, 'd-', label='Total', linewidth=lw, markersize=mks)
    axs[1].set_title('(b) Compositional flow', fontsize=23, pad=10)
    axs[1].set_xlabel('Number of GPUs (Global DOFs)', fontsize=20, labelpad=14)
    axs[1].set_ylabel('Time [s]', fontsize=19)
    axs[1].set_xscale('log', base=2)
    axs[1].set_xticks(ranks2)
    axs[1].set_xticklabels([f"{r:,} ({d})" for r, d in zip(ranks2, global_dofs2)], fontsize=17, rotation=30)
    axs[1].tick_params(axis='y', labelsize=17)
    axs[1].set_ylim(bottom=0)
    axs[1].grid(True, which='both', linestyle='--')

    # Adding a single legend outside the plots
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(1.0, 0.75), fontsize=18)

    # Layout adjustment
    plt.tight_layout(rect=[0, 0, 0.86, 1])  # Adjust the right margin to fit the legend
    print(f"Saving figure {savename}...")
    plt.savefig(savename)
    plt.show()

if __name__ == "__main__":
    weak_scaling_plot("nearwell_scaling_frontier.pdf")
