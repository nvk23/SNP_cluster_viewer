import os
import sys
import subprocess
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.io as pio
import seaborn as sns
from PIL import Image
import datetime
from io import StringIO


def plot_clusters(df, x_col='theta', y_col='r', gtype_col='gt', title='snp plot', opacity = 1, single_snp = None):
    d3 = px.colors.qualitative.D3

    cmap = {
        'AA': d3[0],
        'AB': d3[1],
        'BA': d3[1],
        'BB': d3[2],
        'NC': d3[3]
    }

    # gtypes_list = (df[gtype_col].unique())
    xmin, xmax = df[x_col].min(), df[x_col].max()
    ymin, ymax = df[y_col].min(), df[y_col].max()

    xlim = [xmin-.1, xmax+.1]
    ylim = [ymin-.1, ymax+.1]

    lmap = {'r':'R','theta':'Theta'}
    smap = {'Control':'circle','PD':'diamond-open-dot'}

    if 'R_tightness' in df.columns:
        if single_snp:
            # fig = px.scatter(df.loc[df['snpID'] == single_snp], x=x_col, y=y_col, color=gtype_col, color_discrete_map = cmap, color_continuous_scale=px.colors.sequential.matter, width=650, height=497, labels=lmap, symbol='phenotype', symbol_map=smap, hover_data=['GT', 'R_tightness', 'Theta_tightness']).update_traces(marker_size=100, marker_color="yellow")
            # # fig.add_traces(px.scatter(df.loc[df['snpID'] == single_snp], x=x_col, y=y_col).update_traces(marker_size=100, marker_color="yellow").data)
            # fig.add_traces(px.scatter(df, x=x_col, y=y_col, color=gtype_col, color_discrete_map = cmap, color_continuous_scale=px.colors.sequential.matter, width=650, height=497, labels=lmap, symbol='phenotype', symbol_map=smap, hover_data=['GT', 'R_tightness', 'Theta_tightness']).data)
            # # fig.update_traces(opacity=opacity)   

            # hold_no_snp = df.drop(df[df['snpID'] == single_snp].index)
            fig = fig = px.scatter(df, x=x_col, y=y_col, color=gtype_col, color_discrete_map = cmap, color_continuous_scale=px.colors.sequential.matter, width=650, height=497, labels=lmap, symbol='phenotype', symbol_map=smap, hover_data=['GT', 'R_tightness', 'Theta_tightness'])
            fig.update_traces(opacity = opacity)
            fig.add_traces(px.scatter(df.loc[df['snpID'] == single_snp], x=x_col, y=y_col).update_traces(marker_color="black").data)
        else:
            fig = px.scatter(df, x=x_col, y=y_col, color=gtype_col, color_discrete_map = cmap, color_continuous_scale=px.colors.sequential.matter, width=650, height=497, labels=lmap, symbol='phenotype', symbol_map=smap, hover_data=['GT', 'R_tightness', 'Theta_tightness'])

    else:
        fig = px.scatter(df, x=x_col, y=y_col, color=gtype_col, color_discrete_map=cmap, width=650, height=497, labels=lmap, symbol='phenotype', symbol_map=smap)

    fig.update_xaxes(range=xlim, nticks=10, zeroline=False)
    fig.update_yaxes(range=ylim, nticks=10, zeroline=False)
    
    fig.update_layout(margin=dict(r=76, t=63, b=75))

    # fig.update_layout(legend=dict(orientation="h", yanchor="bottom", y=1, xanchor="right", x=1))

    fig.update_layout(legend_title_text='Genotype')

    out_dict = {
        'fig': fig,
        'xlim': xlim,
        'ylim': ylim
    }
    
    fig.update_layout(title_text=f'<b>{title}<b>')
    
    return out_dict