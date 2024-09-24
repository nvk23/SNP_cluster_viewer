import os
import sys
import subprocess
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.io as pio
import plotly.graph_objects as go
import seaborn as sns
from PIL import Image
import datetime
from io import StringIO


def plot_clusters_seaborn(df, x_col='theta', y_col='r', gtype_col='gt', title='snp plot', opacity=1, highlight_samples=[]):
    # d3 = sns.color_palette("muted")

    # Use a colorblind-friendly palette
    d3 = sns.color_palette("colorblind")
    # d3 = CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
    #               '#f781bf', '#a65628', '#984ea3',
    #               '#999999', '#e41a1c', '#dede00']
    

    # cmap = {
    #     'AA': d3[0],
    #     'AB': d3[1],
    #     'BA': d3[1],
    #     'BB': d3[2],
    #     'NC': d3[3]
    # }

    cmap = {
        'AA': d3[2],
        'AB': d3[4],
        'BA': d3[4],
        'BB': d3[9],
        'NC': d3[7]
    }

    df.dropna(subset=[x_col], inplace = True)
    
    xmin, xmax = df[x_col].min(), df[x_col].max()
    ymin, ymax = df[y_col].min(), df[y_col].max()

    # xlim = [xmin - 0.1, xmax + 0.1]
    xlim = [0, 1.1]
    ylim = [ymin - 0.1, ymax + 0.1]

    lmap = {'r': 'R', 'theta': 'Theta'}
    smap = {'Control': 'o', 'PD': 'D'}

    plt.figure(figsize=(8, 6))

    if len(highlight_samples) > 0:
        # Plot the rest of the data first
        sns.scatterplot(
            data=df[~ df['IID'].isin(highlight_samples)],
            x=x_col,
            y=y_col,
            hue=gtype_col,
            palette=cmap,
            alpha=opacity,
            edgecolor=None,
            legend='full'
        )

        # Plot the single sample on top
        highlight_samples_data = df[df['IID'].isin(highlight_samples)]
        plt.scatter(
            highlight_samples_data[x_col],
            highlight_samples_data[y_col],
            color='black',
            s=200,  # Size of the marker
            edgecolor='black',
            marker='X'  # Diamond D or X marker X
            # label='Selected Samples' # if want in legend
        )
    else:
        # Plot all data
        sns.scatterplot(
            data=df,
            x=x_col,
            y=y_col,
            hue=gtype_col,
            palette=cmap,
            alpha=opacity,
            edgecolor=None,
            legend='full'
        )

    # Set axis limits
    plt.xlim(xlim)
    plt.ylim(ylim)

    # Set axis labels
    plt.xlabel(lmap.get(x_col, x_col))
    plt.ylabel(lmap.get(y_col, y_col))

    # Set title
    plt.title(title)

    # Adjust layout
    plt.legend(title='Genotype', loc='upper left')
    plt.tight_layout()

    # Create the output dictionary
    out_dict = {
        'fig': plt.gcf(),  # Get the current figure
        'xlim': xlim,
        'ylim': ylim
    }
    
    # Adjust the plot size
    plt.gcf().set_size_inches(8, 6)
    plt.gcf().set_figheight(6.25)

    return out_dict