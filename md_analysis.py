import matplotlib.pyplot as plt
import numpy as np
from os.path import exists
import pickle


def plot_md(md_axis, md_data, er_bars, output_fn, show):
    """plots md data with error bars"""
    rd = md_data.shape[0]
    plt.figure(0, figsize=(12, 12))

    # loops over diagonal elements
    for i in range(rd):
        elem = i
        plt.subplot(2, 2, (i + 1))
        plt.errorbar(md_axis, md_data[elem, elem, :], yerr=er_bars[elem, elem, :],
                     ls=" ", mfc="None", marker='o', ecolor='black', color='black',
                     capsize=4, markersize=8, elinewidth=3)
        plt.xlim(0.0, 20.0)
        plt.xticks(fontsize=0.0)
        plt.yticks(fontsize=0.0)

    plt.tight_layout(pad=1.5)

    if not exists(output_fn):
        plt.savefig(output_fn)

    if show is False:
        plt.close()
