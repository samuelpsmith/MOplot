# Import seaborn
# test commits
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import seaborn.objects as so
import os

relative_path = os.path.normpath("\MTPP.csv")
head, tail = os.path.split((os.path.abspath(__file__)))
full_path = head + relative_path
full_path = os.path.normpath(full_path)

#data_path = r'C:\Users\lette\OneDrive\Documents\MTPP.csv'
data_path = full_path
degen = 0.05

def plotMO_swarm(data_path, degen):
    """PlotMO:
    Given a path to a formatted molecular orbital list with labels, plot molecular orbitals. Use Seaborn's swarmplot to reduce overlapping plotting.
    
    """

    ## Set theme and define figure size.
    sns.set_theme(style='ticks', context='paper')
    fig, ax = plt.subplots(figsize=(6.5,12))
    ax.grid(axis='y')
    degen = degen

    ## Define variables
    dataset = pd.read_csv(data_path)
    orbital_num= dataset.loc[:, "orbital_num"]
    compound = dataset.loc[:, "compound"]
    eV = dataset.loc[:, "eV"]
    symmetry_label = dataset.loc[:, "symmetry_label"]
 
    ## Plot data 
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
                #dodge=True
                )

    ## Plot labels and fix labels of degenerate energy levels. 

    for i in range(0, len(compound)):
        if i != 0 and abs(eV[i] - eV[i-1]) <= degen and abs(eV[i] - eV[i-2]) <= degen and abs(eV[i] - eV[i-3]) <= degen and compound[i] == compound[i-1] == compound[i-2] == compound[i-3]:
            # Annotate triply degenerate points with an offset of (100, 0)
            ax.annotate(symmetry_label[i], xy=(compound[i], eV[i]), xytext=(100, 0), size=10,
                        ha="center", va="top", textcoords="offset points")
        elif i != 0 and abs(eV[i] - eV[i-1]) <= degen and abs(eV[i] - eV[i-2]) <= degen and compound[i] == compound[i-1] == compound[i-2]:
            # Annotate triply degenerate points with an offset of (80, 0)
            ax.annotate(symmetry_label[i], xy=(compound[i], eV[i]), xytext=(80, 0), size=10,
                        ha="center", va="top", textcoords="offset points")
        elif i != 0 and abs(eV[i] - eV[i-1]) <= degen and compound[i] == compound[i-1]:
            # Annotate doubly degenerate points with an offset of (60, 0)
            ax.annotate(symmetry_label[i], xy=(compound[i], eV[i]), xytext=(60, 0), size=10,
                        ha="center", va="top", textcoords="offset points")
        else:
            # Annotate non-degenerate points with an offset of (40, 0)
            ax.annotate(symmetry_label[i], xy=(compound[i], eV[i]), xytext=(40, 0), size=10,
                        ha="center", va="top", textcoords="offset points")

    sns.despine()
    plt.show()

def plotMO_cat(data_path, degen):
    """PlotMO:
    Given a path to a formatted molecular orbital list with labels, plot molecular orbitals.
    
    """
    
    ## Set theme and define figure size.
    sns.set_theme(style='ticks', context='paper')
    fig, ax = plt.subplots(figsize=(6.5,12))
    ax.grid(axis='y')
    degen = degen

    ## Define variables
    dataset = pd.read_csv(data_path)
    orbital_num= dataset.loc[:, "orbital_num"]
    compound = dataset.loc[:, "compound"]
    eV = dataset.loc[:, "eV"]
    symmetry_label = dataset.loc[:, "symmetry_label"]
 
    ## Plot data 
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
                #dodge=True
                )

    ## Plot labels and fix labels of degenerate energy levels. 

    for i in range(0, len(compound)):
        if i != 0 and abs(eV[i] - eV[i-1]) <= degen and abs(eV[i] - eV[i-2]) <= degen and abs(eV[i] - eV[i-3]) <= degen and compound[i] == compound[i-1] == compound[i-2] == compound[i-3]:
            # Annotate triply degenerate points with an offset of (100, 0)
            ax.annotate(symmetry_label[i], xy=(compound[i], eV[i]), xytext=(100, 0), size=10,
                        ha="center", va="top", textcoords="offset points")
        elif i != 0 and abs(eV[i] - eV[i-1]) <= degen and abs(eV[i] - eV[i-2]) <= degen and compound[i] == compound[i-1] == compound[i-2]:
            # Annotate triply degenerate points with an offset of (80, 0)
            ax.annotate(symmetry_label[i], xy=(compound[i], eV[i]), xytext=(80, 0), size=10,
                        ha="center", va="top", textcoords="offset points")
        elif i != 0 and abs(eV[i] - eV[i-1]) <= degen and compound[i] == compound[i-1]:
            # Annotate doubly degenerate points with an offset of (60, 0)
            ax.annotate(symmetry_label[i], xy=(compound[i], eV[i]), xytext=(60, 0), size=10,
                        ha="center", va="top", textcoords="offset points")
        else:
            # Annotate non-degenerate points with an offset of (40, 0)
            ax.annotate(symmetry_label[i], xy=(compound[i], eV[i]), xytext=(40, 0), size=10,
                        ha="center", va="top", textcoords="offset points")

    sns.despine()
    plt.show()

plotMO_swarm(data_path, degen)