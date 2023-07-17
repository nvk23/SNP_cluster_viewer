import os
import sys
import subprocess
import datetime
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.io as pio
import seaborn as sns
import random
from PIL import Image
from st_aggrid import GridOptionsBuilder, AgGrid, GridUpdateMode, DataReturnMode
from io import StringIO, BytesIO

from hold_method import plot_clusters

st.title('Interactive Previous No Calls Evaluation')

if 'snp_choices' not in st.session_state:
    st.session_state['snp_choices'] = []

full_metrics = pd.read_csv('data/060623_full_cluster_tightness_wip')
tightness_measure = st.sidebar.selectbox(label = 'Choose the Outlier Metric', options=['R', 'Theta'])
rand_generate = st.sidebar.checkbox('Randomly generate SNPs')

if tightness_measure == 'R':
    measure = 'R_outlier'
elif tightness_measure == 'Theta':
    measure = 'Theta_outlier'

flagged_outliers = full_metrics[full_metrics[measure] == True]
prev_nc = flagged_outliers[flagged_outliers["GT"] == "NC"]
now_nc = flagged_outliers[flagged_outliers["preds_cat"] == "NC"]

exclude_snps = set(now_nc['snpID'].unique())
include_snps = set(prev_nc['snpID'].unique())
include = include_snps.difference(exclude_snps)

if not rand_generate:
    st.sidebar.markdown('### Choose an outlier SNP to display')
    snp_name = st.sidebar.selectbox(label = 'SNP Name Choice', label_visibility = 'collapsed', options=include)

    snp1 = full_metrics.loc[full_metrics['snpID'] == snp_name]
else:
    snp_name = random.choice(list(include))
    snp1 = full_metrics.loc[full_metrics['snpID'] == snp_name]

    st.markdown(f'#### SNP: {snp_name}')
    btn1, btn2, btn3 = st.columns([1, 0.5, 0.5])
    btn1.markdown('#### Should this reclustering be regenerated?')
    regenerate = btn2.button('Yes')
    no_regenerate = btn3.button('No')

    if regenerate:
        st.session_state['snp_choices'].append(snp_name)
        snp_name = random.choice(list(include))
        snp1 = full_metrics.loc[full_metrics['snpID'] == snp_name]
    elif no_regenerate:
        st.session_state['snp_choices'].append(snp_name)
        snp_name = random.choice(list(include))
        snp1 = full_metrics.loc[full_metrics['snpID'] == snp_name] 

    st.sidebar.markdown('Reported SNPs:')
    st.sidebar.write(st.session_state['snp_choices'])  # MAKE DOWNLOAD BUTTON FOR TXT FILE OF THESE SNPS

cluster_metric_before = plot_clusters(snp1, x_col='Theta', y_col='R', gtype_col='GT', title = "Before Recluster")['fig']
st.plotly_chart(cluster_metric_before, use_container_width = True)

cluster_metric_after = plot_clusters(snp1, x_col='Theta', y_col='R', gtype_col='preds_cat', title = "After Recluster")['fig']
st.plotly_chart(cluster_metric_after, use_container_width=True)

cluster_full = plot_clusters(full_metrics, x_col='Theta', y_col='R', gtype_col='preds_cat', title = "Full Plot", single_snp = snp_name, opacity = 0.2)['fig']
st.plotly_chart(cluster_full, use_container_width=True)

