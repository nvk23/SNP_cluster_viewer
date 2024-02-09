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
from PIL import Image
from st_aggrid import GridOptionsBuilder, AgGrid, GridUpdateMode, DataReturnMode
from io import StringIO, BytesIO

from hold_method import plot_clusters

st.title('Previous No Calls Evaluation')

st.sidebar.markdown('### Choose a Model')
models = ['model_060623', 'model_080823_allgentrainscores', 'model_081623_small', 'model_081623']
model_name = st.sidebar.selectbox(label = 'Model Choice', label_visibility = 'collapsed', options=models)

full_metrics = pd.read_csv(f'data/{model_name}_full_cluster_tightness.csv')

tightness_measure = st.sidebar.selectbox(label = 'Tightness Metric', options=['None', 'R', 'Theta'])

if tightness_measure == 'R':
    r_less = full_metrics.loc[full_metrics['R_tightness']<1.042037]
    prev_nc = r_less[r_less["GT"] == "NC"]
    now_nc = r_less[r_less["preds_cat"] == "NC"]
    # st.dataframe(now_nc)

    exclude_snps = set(now_nc['snpID'].unique())
    include_snps = set(prev_nc['snpID'].unique())
    include = include_snps.difference(exclude_snps)

    st.sidebar.markdown('### Choose an individual SNP to display')
    snp_name = st.sidebar.selectbox(label = 'SNP Name Choice', label_visibility = 'collapsed', options=include)

    snp1 = full_metrics.loc[full_metrics['snpID'] == snp_name]

    cluster_metric_before = plot_clusters(snp1, x_col='Theta', y_col='R', gtype_col='GT', title = "Before Recluster")['fig']
    st.plotly_chart(cluster_metric_before, use_container_width = True)

    cluster_metric_after = plot_clusters(snp1, x_col='Theta', y_col='R', gtype_col='preds_cat', title = "After Recluster")['fig']
    st.plotly_chart(cluster_metric_after, use_container_width=True)

    cluster_full = plot_clusters(full_metrics, x_col='Theta', y_col='R', gtype_col='preds_cat', title = "After Recluster", single_snp = snp_name, opacity = 0.2)['fig']
    st.plotly_chart(cluster_full, use_container_width=True)

elif tightness_measure == 'Theta':
    theta_less = full_metrics.loc[full_metrics['Theta_tightness']<0.059508]

    prev_nc = theta_less[theta_less["GT"] == "NC"]
    now_nc = theta_less[theta_less["preds_cat"] == "NC"]
    # st.dataframe(now_nc)

    exclude_snps = set(now_nc['snpID'].unique())
    include_snps = set(prev_nc['snpID'].unique())
    include = include_snps.difference(exclude_snps)

    st.sidebar.markdown('### Choose an individual SNP to display')
    snp_name = st.sidebar.selectbox(label = 'SNP Name Choice', label_visibility = 'collapsed', options=include)

    snp1 = full_metrics.loc[full_metrics['snpID'] == snp_name]

    cluster_metric_before = plot_clusters(snp1, x_col='Theta', y_col='R', gtype_col='GT', title = "Before Recluster")['fig']
    st.plotly_chart(cluster_metric_before, use_container_width = True)

    cluster_metric_after = plot_clusters(snp1, x_col='Theta', y_col='R', gtype_col='preds_cat', title = "After Recluster")['fig']
    st.plotly_chart(cluster_metric_after, use_container_width=True)

    cluster_full = plot_clusters(full_metrics, x_col='Theta', y_col='R', gtype_col='preds_cat', title = "After Recluster", single_snp = snp_name, opacity = 0.2)['fig']
    st.plotly_chart(cluster_full, use_container_width=True)


elif tightness_measure == 'None':
    prev_nc = full_metrics[full_metrics["GT"] == "NC"]
    now_nc = full_metrics[full_metrics["preds_cat"] == "NC"]
    # st.dataframe(now_nc)

    exclude_snps = set(now_nc['snpID'].unique())
    include_snps = set(prev_nc['snpID'].unique())
    include = include_snps.difference(exclude_snps)

    st.sidebar.markdown('### Choose an individual SNP to display')
    snp_name = st.sidebar.selectbox(label = 'SNP Name Choice', label_visibility = 'collapsed', options=include)

    snp1 = full_metrics.loc[full_metrics['snpID'] == snp_name]

    cluster_metric_before = plot_clusters(snp1, x_col='Theta', y_col='R', gtype_col='GT', title = "Before Recluster")['fig']
    st.plotly_chart(cluster_metric_before, use_container_width = True)

    cluster_metric_after = plot_clusters(snp1, x_col='Theta', y_col='R', gtype_col='preds_cat', title = "After Recluster")['fig']
    st.plotly_chart(cluster_metric_after, use_container_width=True)

    cluster_full = plot_clusters(full_metrics, x_col='Theta', y_col='R', gtype_col='preds_cat', title = "After Recluster", single_snp = snp_name, opacity = 0.2)['fig']
    st.plotly_chart(cluster_full, use_container_width=True)

