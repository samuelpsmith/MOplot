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
    relative_path = os.path.normpath("/MTPFPP.csv")
    head, tail = os.path.split((os.path.abspath(__file__)))
    full_path = head + relative_path
    full_path = os.path.normpath(full_path)
    return full_path


def plotMO_swarm(data_path, degen, size, figsize, textX, marker_size, line_width, texty):
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

    #Apply a small vertical offset to the data points to avoid overlapping points
    #dataset['eV'] = dataset['eV'] + np.random.uniform(-0.1, 0.1, len(dataset))

    # Plot data
    sns.swarmplot(data=dataset, x='compound', y='eV',
                  marker='_',
                  linewidth=line_width, # width of markers
                  zorder=3,
                  s=marker_size, # size of markers
                  hue='symmetry_label',
                  edgecolor='face',
                  palette='flare_r',
                  ax=ax,
                  legend=False,
                  dodge=False
                  )
    ax.set(xlabel=None)

    # Plot labels and fix labels of degenerate energy levels

    # Define offsets for labels, based on textX:
    textX = float(textX)
    # produce the variables offD3, offD2, offD, offset without explicitly casting them as floats:
    #offD3 = (8*textX, 0)
    #offD2 = (6*textX, 0)
    #offD = (4*textX, 0)
    #offset = (2*textX, 0)
    
    offD3 = (0.5 * marker_size + 8*textX, line_width + texty)
    offD2 = (0.5 * marker_size + 6*textX, line_width + texty)
    offD = (0.5 * marker_size + 4*textX, line_width + texty)
    offset = (0.5 * marker_size + 2*textX, line_width + texty)

    degen_num = check_degen(eV, degen)


    for i in range(0, len(compound)):
        if degen_num[i] == 0:
            # Annotate non-degenerate points with an offset of offset
            ax.annotate(symmetry_label[i], xy=(compound[i], eV[i]), xytext=offset, size=size,
                        ha="center", va="top", textcoords="offset points")
        elif degen_num[i] == 1:
            # Annotate degenerate points with an offset of offD
            ax.annotate(symmetry_label[i], xy=(compound[i], eV[i]), xytext=offD, size=size,
                        ha="center", va="top", textcoords="offset points")
        elif degen_num[i] == 2:
            # Annotate doubly degenerate points with an offset of offD2
            ax.annotate(symmetry_label[i], xy=(compound[i], eV[i]), xytext=offD2, size=size,
                        ha="center", va="top", textcoords="offset points")
        elif degen_num[i] == 3:
            # Annotate triply degenerate points with an offset of offD3
            ax.annotate(symmetry_label[i], xy=(compound[i], eV[i]), xytext=offD3, size=size,
                        ha="center", va="top", textcoords="offset points")
        else:
            print("Error: check_degen() returned a value greater than 3.")


    sns.despine()
    plt.show()

# for each entry in an array, check how many previous entries are within some tolerance of it


def check_degen(array, degen):
    """check_degen:
    Given an array of values, check how many previous entries are within some tolerance of it and precede it by less than 4 places. Return an array of integers with the same length as the input array.
    """
    # consider changing this so that only things with the same x axis coordinate are considered degenerate
    
    degen = degen
    degen_array = np.zeros(len(array), dtype=int)
    for i in range(0, len(array)):
        for j in range(0, i):
            if (abs(array[i] - array[j]) < degen) and i - j < 4:
                degen_array[i] += 1
    return degen_array


def plotMO_cat(data_path, degen, size, figsize, textX, marker_size, line_width, texty, vertical_jitter):
    """PlotMO:
    Given a path to a formatted molecular orbital list with labels, plot molecular orbitals.

    """

    # Set theme and define figure size.
    sns.set_theme(style='ticks', context='paper')
    fig, ax = plt.subplots(figsize=figsize)
    ax.grid(axis='y')
    #degen = degen

    # Define variables
    dataset = pd.read_csv(data_path)
    #orbital_num = dataset.loc[:, "orbital_num"]
    compound = dataset.loc[:, "compound"]
    eV = dataset.loc[:, "eV"]
    symmetry_label = dataset.loc[:, "symmetry_label"]

    #Apply a small vertical offset to the data points to avoid overlapping points
    np.random.seed(2)
    degen_num = check_degen(eV, degen)
    for i in range(0, len(compound)):
        if degen_num[i] > 0:
            # Apply small vetical offset to degenerate points\
            dataset['eV'][i] = dataset['eV'][i] + np.random.uniform(-vertical_jitter, vertical_jitter, 1)
            
    #dataset['eV'] = dataset['eV'] + np.random.uniform(-0.1, 0.1, len(dataset))
      

    # Plot data
    sns.stripplot(data=dataset, x='compound', y='eV',
                  marker='_',
                  linewidth=line_width, # width of markers
                  zorder=3,
                  s=marker_size, # size of markers
                  hue='symmetry_label',
                  edgecolor='face',
                  palette='flare_r',
                  ax=ax,
                  legend=False,
                  jitter=False,
                  dodge=False
                  )
    

    ax.set(xlabel=None)

    # Plot labels and fix labels of degenerate energy levels

    # Define offsets for labels, based on textX:
    textX = float(textX)
    # produce the variables offD3, offD2, offD, offset without explicitly casting them as floats:
    #offD3 = (8*textX, 0)
    #offD2 = (6*textX, 0)
    #offD = (4*textX, 0)
    #offset = (2*textX, 0)
    
    offD3 = (0.5 * marker_size + 8*textX, line_width + texty)
    offD2 = (0.5 * marker_size + 6*textX, line_width + texty)
    offD = (0.5 * marker_size + 4*textX, line_width + texty)
    offset = (0.5 * marker_size + 2*textX, line_width + texty)

    degen_num = check_degen(eV, degen)

    for i in range(0, len(compound)):
        if degen_num[i] == 0:
            # Annotate non-degenerate points with an offset of offset
            ax.annotate(symmetry_label[i], xy=(compound[i], eV[i]), xytext=offset, size=size,
                        ha="center", va="top", textcoords="offset points")
        elif degen_num[i] == 1:
            # Annotate degenerate points with an offset of offD
            ax.annotate(symmetry_label[i], xy=(compound[i], eV[i]), xytext=offD, size=size,
                        ha="center", va="top", textcoords="offset points")
        elif degen_num[i] == 2:
            # Annotate doubly degenerate points with an offset of offD2
            ax.annotate(symmetry_label[i], xy=(compound[i], eV[i]), xytext=offD2, size=size,
                        ha="center", va="top", textcoords="offset points")
        elif degen_num[i] == 3:
            # Annotate triply degenerate points with an offset of offD3
            ax.annotate(symmetry_label[i], xy=(compound[i], eV[i]), xytext=offD3, size=size,
                        ha="center", va="top", textcoords="offset points")
        else:
            print("Error: check_degen() returned a value greater than 3.")

    # change axis label size using plt.xlabel
    #plt.xlabel('Compound', fontsize=14)
    plt.ylabel('eV', fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    sns.despine()
    plt.show()


data_path = get_path() # path to data file
degen = 0.05 # degeneracy tolerance
size = 11 # size of labels
figsize = (6.5, 4.5) # size of figure
textX = 11 # multiple for x offset of labels
texty = 2 #y offset of labels added to marker hight
marker_size = 20 # size of markers
line_width = 3 # width of markers
vertical_jitter = 0.075 # amount of vertical jitter applied to degenerate points
plotMO_cat(data_path, degen, size, figsize, textX, marker_size, line_width, texty, vertical_jitter)
#plotMO_swarm(data_path, degen, size, figsize, textX, marker_size, line_width, texty)
