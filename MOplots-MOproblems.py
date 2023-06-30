# To Do:
# - Add option to plot with or without labels

# Import modules
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os


def get_path():
    """get_path:
    Return the full path to the data file, regardless of where the script is run from.
    """
    relative_path = os.path.normpath("/MTPP.csv")
    head, tail = os.path.split((os.path.abspath(__file__)))
    full_path = head + relative_path
    full_path = os.path.normpath(full_path)
    return full_path


def plotMO_swarm(data_path, degen, size, figsize, textX):
    """PlotMO:
    Given a path to a formatted molecular orbital list with labels, plot molecular orbitals. Use Seaborn's swarmplot to reduce overlapping plotting.

    """

    # Set theme and define figure size.
    sns.set_theme(style='ticks', context='paper')
    fig, ax = plt.subplots(figsize=figsize)
    ax.grid(axis='y')
    degen = degen

    # Define variables
    dataset = pd.read_csv(data_path)
    orbital_num = dataset.loc[:, "orbital_num"]
    compound = dataset.loc[:, "compound"]
    eV = dataset.loc[:, "eV"]
    symmetry_label = dataset.loc[:, "symmetry_label"]

    # Plot data
    sns.swarmplot(data=dataset, x='compound', y='eV',
                  marker='_',
                  linewidth=1,
                  zorder=3,
                  s=20,
                  hue='orbital_label',
                  edgecolor='face',
                  palette='rocket',
                  ax=ax,
                  legend=False,
                  # dodge=True
                  )
    ax.set(xlabel=None)

    # for each entry in an array, check how many previous entries are within degen of it
    for i in range(0, len(compound)):
        degen_num = check_degen(eV, degen)
        if degen_num[i] == 0:
            # Annotate non-degenerate points with an offset of (40, 0)
            ax.annotate(symmetry_label[i], xy=(compound[i], eV[i]), xytext=(40, 0), size=size,
                        ha="center", va="top", textcoords="offset points")
        elif degen_num[i] == 1:
            # Annotate doubly degenerate points with an offset of (60, 0)
            ax.annotate(symmetry_label[i], xy=(compound[i], eV[i]), xytext=(60, 0), size=size,
                        ha="center", va="top", textcoords="offset points")
        elif degen_num[i] == 2:
            # Annotate triply degenerate points with an offset of (80, 0)
            ax.annotate(symmetry_label[i], xy=(compound[i], eV[i]), xytext=(80, 0), size=size,
                        ha="center", va="top", textcoords="offset points")
        elif degen_num[i] == 3:
            # Annotate triply degenerate points with an offset of (100, 0)
            ax.annotate(symmetry_label[i], xy=(compound[i], eV[i]), xytext=(100, 0), size=size,
                        ha="center", va="top", textcoords="offset points")

    sns.despine()
    plt.show()

# for each entry in an array, check how many previous entries are within some tolerance of it


def check_degen(array, degen):
    """check_degen:
    Given an array of values, check how many previous entries are within some tolerance of it. Return an array of integers with the same length as the input array.
    """
    degen = degen
    degen_array = np.zeros(len(array), dtype=int)
    for i in range(0, len(array)):
        for j in range(0, i):
            if abs(array[i] - array[j]) < degen:
                degen_array[i] += 1
    return degen_array


def plotMO_cat(data_path, degen, size, figsize, textX):
    """PlotMO:
    Given a path to a formatted molecular orbital list with labels, plot molecular orbitals.

    """

    # Set theme and define figure size.
    sns.set_theme(style='ticks', context='paper')
    fig, ax = plt.subplots(figsize=figsize)
    ax.grid(axis='y')
    degen = degen

    # Define variables
    dataset = pd.read_csv(data_path)
    orbital_num = dataset.loc[:, "orbital_num"]
    compound = dataset.loc[:, "compound"]
    eV = dataset.loc[:, "eV"]
    symmetry_label = dataset.loc[:, "symmetry_label"]

    # Plot data
    sns.stripplot(data=dataset, x='compound', y='eV',
                  marker='_',
                  linewidth=1,
                  zorder=3,
                  s=20,
                  hue='orbital_label',
                  edgecolor='face',
                  palette='rocket',
                  ax=ax,
                  legend=False,
                  jitter=False
                  # dodge=True
                  )
    ax.set(xlabel=None)

    # Plot labels and fix labels of degenerate energy levels

    # Define offsets for labels, based on textX:
    textX = float(textX)
    # produce the variables offD3, offD2, offD, offset without explicitly casting them as floats:
    offD3 = (10*textX, 0)
    offD2 = (8*textX, 0)
    offD = (6*textX, 0)
    offset = (4*textX, 0)

    print(len(compound))

    for i in range(0, len(compound)):
        degen_num = check_degen(eV, degen)
        if degen_num[i] == 0:
            # Annotate non-degenerate points with an offset of offset
            ax.annotate(symmetry_label[i], xy=(compound[i], eV[i]), xytext=offset, size=size,
                        ha="center", va="top", textcoords="offset points")
        elif degen_num[i] == 1:
            # Annotate doubly degenerate points with an offset of (60, 0)
            ax.annotate(symmetry_label[i], xy=(compound[i], eV[i]), xytext=offD, size=size,
                        ha="center", va="top", textcoords="offset points")
        elif degen_num[i] == 2:
            # Annotate triply degenerate points with an offset of (80, 0)
            ax.annotate(symmetry_label[i], xy=(compound[i], eV[i]), xytext=offD2, size=size,
                        ha="center", va="top", textcoords="offset points")
        elif degen_num[i] == 3:
            # Annotate triply degenerate points with an offset of (100, 0)
            ax.annotate(symmetry_label[i], xy=(compound[i], eV[i]), xytext=offD3, size=size,
                        ha="center", va="top", textcoords="offset points")
        else:
            print("Error: check_degen() returned a value greater than 3.")

    sns.despine()
    plt.show()


data_path = get_path()
degen = 0.03
size = 10
figsize = (6, 4)
textX = 10
plotMO_cat(data_path, degen, size, figsize, textX)
