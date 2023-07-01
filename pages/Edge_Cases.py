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
from Metrics_Home import plot_clusters

cluster_metrics = pd.read_csv(f'data/model_060623_full_cluster_tightness')
tightness_measure = st.sidebar.selectbox(label = 'Cluster Tightnes Metric', options=['R', 'Theta'])

percentile = st.sidebar.checkbox('Show percentiles', value = True)
metric_level = st.sidebar.checkbox('Choose metric level')
min_r = 0.174188
max_r = 62.392498
min_theta = 0.007535
max_theta = 48.659265

if metric_level:
    if tightness_measure == 'R':
        metric_value = st.select_slider('Display R tightness under this value:', options=np.arange(min_r, max_r, 0.5).tolist())
    elif tightness_measure == 'Theta':
        metric_value = st.select_slider('Display Theta tightness under this value:', options=np.arange(min_theta, max_theta, 0.5).tolist())

    col_name = f'{tightness_measure}_tightness'
    title = f'{tightness_measure} Tightness Under {metric_value}'
    cluster_custom = cluster_metrics.loc[cluster_metrics[col_name]<metric_value]
    cluster_metric = plot_clusters(cluster_custom, x_col='Theta', y_col='R', gtype_col='preds_cat', title = title)['fig']
    st.plotly_chart(cluster_metric, use_container_width=True)

    title = f'{tightness_measure} Tightness Over {metric_value}'
    cluster_custom = cluster_metrics.loc[cluster_metrics[col_name]>metric_value]
    cluster_metric = plot_clusters(cluster_custom, x_col='Theta', y_col='R', gtype_col='preds_cat', title = title)['fig']
    st.plotly_chart(cluster_metric, use_container_width=True)

if percentile: 
    cluster_25_r = cluster_metrics.loc[cluster_metrics['R_tightness']<1.042037]
    cluster_75_r = cluster_metrics.loc[cluster_metrics['R_tightness']>3.122486]

    cluster_25_theta = cluster_metrics.loc[cluster_metrics['Theta_tightness']<0.059508]
    cluster_75_theta = cluster_metrics.loc[cluster_metrics['Theta_tightness']>1.156379]

    if tightness_measure == 'R':
        cluster_25 = plot_clusters(cluster_25_r, x_col='Theta', y_col='R', gtype_col='preds_cat', title = '25th Percentile R Tightness')['fig']
        st.plotly_chart(cluster_25, use_container_width=True)

        cluster_25 = plot_clusters(cluster_75_r, x_col='Theta', y_col='R', gtype_col='preds_cat', title = 'Over 75th Percentile R Tightness')['fig']
        st.plotly_chart(cluster_25, use_container_width=True)

    elif tightness_measure == 'Theta':
        cluster_25 = plot_clusters(cluster_25_theta, x_col='Theta', y_col='R', gtype_col='preds_cat', title = '25th Percentile Theta Tightness')['fig']
        st.plotly_chart(cluster_25, use_container_width=True)

        cluster_25 = plot_clusters(cluster_75_theta, x_col='Theta', y_col='R', gtype_col='preds_cat', title = 'Over 75th Percentile Theta Tightness')['fig']
        st.plotly_chart(cluster_25, use_container_width=True)