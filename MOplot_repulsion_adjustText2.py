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

from adjustText import adjust_text

def repel_annotations(ax, annotations):
    """
    Detect and resolve overlapping annotations on a given Axes object using adjustText.

    :param ax: The Axes object containing the annotations.
    :param annotations: List of the annotations on the Axes.
    :return: None
    """
    # List to store texts for adjust_text
    texts = []
    
    # Extract the text objects from the annotations
    for ann in annotations:
        # Create a text object with the same properties as the annotation
        t = ax.text(ann.get_position()[0], ann.get_position()[1], ann.get_text(),
                    ha=ann.get_horizontalalignment(),
                    va=ann.get_verticalalignment(),
                    fontsize=ann.get_fontsize(),
                    color=ann.get_color())
        texts.append(t)
        
        # Remove the original annotations as we're now using text objects
        ann.remove()

    # Use adjust_text to repel overlapping texts
    adjust_text(texts)

    # Redraw the figure after repelling texts
    ax.figure.canvas.draw_idle()

def points_to_data_x(ax, points_distance_x):
    # Transform a point at the origin to display coordinates
    disp_coord = ax.transData.transform((0, 0))
    
    # Add the desired offset in points
    disp_coord_offset = disp_coord[0] + points_distance_x, disp_coord[1]
    
    # Inverse transform the new display coordinate back to data coordinates
    data_coord = ax.transData.inverted().transform(disp_coord_offset)
    
    # Return the difference, which represents the points_distance in data coordinates
    return data_coord[0] - 0

def points_to_data_y(ax, points_distance_y):
    # Transform a point at the origin to display coordinates
    disp_coord = ax.transData.transform((0, 0))
    
    # Add the desired offset in points
    disp_coord_offset = disp_coord[0], disp_coord[1] + points_distance_y
    
    # Inverse transform the new display coordinate back to data coordinates
    data_coord = ax.transData.inverted().transform(disp_coord_offset)
    
    # Return the difference, which represents the points_distance in data coordinates
    return data_coord[1] - 0

def points_to_data(ax, offset):
    return points_to_data_x(ax, offset[0]), points_to_data_y(ax, offset[1])


def plotMO_cat(data_path, degen, size, figsize, textX, marker_size, line_width, texty, vertical_jitter, use_adjustText):
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
        if degen_num[i] == 0:
            # Apply no offset for nondegenerate points
            pass
        elif degen_num[i] == 1:
            # Annotate degenerate points with an offset of offD
            dataset['eV'][i] = dataset['eV'][i] + np.random.uniform(-vertical_jitter, vertical_jitter, 1)
        elif degen_num[i] == 2:
            # Annotate doubly degenerate points with an offset of offD2
            dataset['eV'][i] = dataset['eV'][i] + 1.05*np.random.uniform(-vertical_jitter, vertical_jitter, 1)
        elif degen_num[i] == 3:
            # Annotate triply degenerate points with an offset of offD3
            dataset['eV'][i] = dataset['eV'][i] + 1.10*np.random.uniform(-vertical_jitter, vertical_jitter, 1)
        else:
            print("Error: check_degen() returned a value greater than 3.")
            
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
    offD3 = (0.5 * marker_size + 8*textX, line_width + texty)
    offD2 = (0.5 * marker_size + 6*textX, line_width + texty)
    offD = (0.5 * marker_size + 4*textX, line_width + texty)
    offset = (0.5 * marker_size + 2*textX, line_width + texty)

    offD3_data = points_to_data(ax, offD3)
    offD2_data = points_to_data(ax, offD2)
    offD_data = points_to_data(ax, offD)
    offset_data = points_to_data(ax, offset)
    

    degen_num = check_degen(eV, degen)

    # Create a mapping of compound names to their numeric x-values
    unique_compounds = dataset['compound'].unique()
    # Recreate the compound_to_x mapping
    compound_to_x = {compound: i for i, compound in enumerate(unique_compounds)}
    print(compound_to_x)
    # Use the mapping to get numeric x-values for each compound in the dataset
    compound_numeric = [compound_to_x[comp] for comp in compound]

    ax.set_xlim(-0.5, len(unique_compounds))


    texts = []  # List to store the text objects

    for i in range(0, len(compound)):
        x = compound_to_x[compound[i]]  # retrieve the numeric x position
        print(x)
        y = eV[i]  # original y position
        label = symmetry_label[i]  # the label text
        print(label)
        #offset = None  # placeholder for the chosen offset
        
        # Assign offsets based on degeneracy
        if degen_num[i] == 0:
            offset = offset_data
        elif degen_num[i] == 1:
            offset = offD_data
        elif degen_num[i] == 2:
            offset = offD2_data
        elif degen_num[i] == 3:
            offset = offD3_data
        else:
            print("Error: check_degen() returned a value greater than 3.")
        
        # Calculate the new x and y positions using the offsets
        new_x = x + offset[0]
        new_y = y + offset[1]
        print(new_x, new_y)
        
        # Use ax.text to add the text at the new position
        text = ax.text(new_x, new_y, label, size=size, ha="center", va="top")

        if use_adjust_text:
            texts.append(text)

    if use_adjust_text: # Apply adjust_text() outside the loop
        adjust_text(texts, x=[t.get_position()[0] for t in texts], y=[t.get_position()[1] for t in texts], 
                    force_text=(0.1, 0.2),
                    force_static=(0.1, 0.2),
                    force_pull=(0.05, 0.01),
                    force_explode=(0.01, 0.02),
                    expand=(1.05, 1.1),
                    ensure_inside_axes=True,
                    ax=ax,
                    avoid_self=False,
                    only_move={'static':'x', 'text':'y'},
                    expand_axes=True)


    # Call repel_annotations
    #annotations = [child for child in ax.get_children() if isinstance(child, plt.Text)]
    #repel_annotations(ax, annotations)
    
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
textX = 15 # multiple for x offset of labels
texty = 2 #y offset of labels added to marker hight
marker_size = 20 # size of markers
line_width = 3 # width of markers
vertical_jitter = 0.025 #5 # amount of vertical jitter applied to degenerate points
use_adjust_text = False # use adjust_text() to avoid overlapping labels
plotMO_cat(data_path, degen, size, figsize, textX, marker_size, line_width, texty, vertical_jitter, use_adjust_text)
#plotMO_swarm(data_path, degen, size, figsize, textX, marker_size, line_width, texty)
