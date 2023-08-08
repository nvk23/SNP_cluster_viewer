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

def find_unmatched(df_all):
    df_all = df_all.reset_index()
    mask = (df_all['GT'] != 'NC') & (df_all['preds_cat'] != 'NC')

    # Check if 'GT' is equal to 'preds' where the mask is True
    result = (df_all['GT'] == df_all['preds_cat'])[mask]
    df_all['not_NC_check'] = result
    unmatched = df_all.loc[df_all['not_NC_check'] == False, 'snpID']
    # display(df_all[['GT', 'R', 'preds_cat', 'not_NC_check']].loc[df_all['not_NC_check'] == False])

    return unmatched.values

def display_plots(snp_name, full_metrics, full_plot = True, before_after = True, title = None): 
    snp1 = full_metrics.loc[full_metrics['snpID'] == snp_name]

    if before_after:
        cluster_metric_before = plot_clusters(snp1, x_col='Theta', y_col='R', gtype_col='GT', title = "Before Recluster")['fig']
        st.plotly_chart(cluster_metric_before, use_container_width = True)

        cluster_metric_after = plot_clusters(snp1, x_col='Theta', y_col='R', gtype_col='preds_cat', title = "After Recluster")['fig']
        st.plotly_chart(cluster_metric_after, use_container_width=True)
    else:
        cluster_metric = plot_clusters(snp1, x_col='Theta', y_col='R', gtype_col='preds_cat', title = title)['fig']
        st.plotly_chart(cluster_metric, use_container_width = True)

    if full_plot:
        cluster_full = plot_clusters(full_metrics, x_col='Theta', y_col='R', gtype_col='preds_cat', title = "Full Plot", single_snp = snp_name, opacity = 0.2)['fig']
        st.plotly_chart(cluster_full, use_container_width=True)


if 'snp_choices' not in st.session_state:
    st.session_state['snp_choices'] = []
if 'snp_seen' not in st.session_state:
    st.session_state['snp_seen'] = []
if 'snp_name' not in st.session_state:
    st.session_state['snp_name'] = ''

full_metrics = pd.read_csv('data/060623_full_cluster_tightness_wip')
before_matched = pd.read_csv('data/060623_complete')

display_unmatched = st.sidebar.checkbox('Display discrepancy SNPs')

if not display_unmatched:
    rand_generate = st.sidebar.checkbox('Randomly generate SNPs')
    tightness_measure = st.sidebar.selectbox(label = 'Choose the Outlier Metric:', options=['R', 'Theta'])

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

        display_plots(snp_name, full_metrics)

    else:
        btn1, btn2, btn3 = st.columns([1, 0.5, 0.5])
        btn1.markdown('#### Should this reclustering be regenerated?')
        regenerate = btn2.button('Yes')
        no_regenerate = btn3.button('No')

        if regenerate:
            st.session_state['snp_choices'].append(st.session_state['snp_name'])
            st.session_state['snp_seen'].append(st.session_state['snp_name'])
            st.session_state['snp_name'] = random.choice(list(include))
        elif no_regenerate:
            st.session_state['snp_seen'].append(st.session_state['snp_name'])
            st.session_state['snp_name'] = random.choice(list(include))
        else:
            st.session_state['snp_name'] = random.choice(list(include))

        try:
            if st.session_state['snp_name'] in st.session_state['snp_seen']:
                st.session_state['snp_name'] = random.choice(list(include))
        except:
            st.error("No more SNPs to run through!")

        st.markdown(f'#### SNP: {st.session_state["snp_name"]}')
        st.sidebar.markdown('Reported SNPs:')
        
        # MAKE DOWNLOAD BUTTON FOR TXT FILE OF THESE SNPS
        st.sidebar.data_editor(st.session_state['snp_choices'],
                        column_config={"value": st.column_config.TextColumn("Outlier SNP")},
                        hide_index=True,
                        use_container_width=True
                    )
        display_plots(st.session_state['snp_name'], full_metrics)

elif display_unmatched:
    # Temporarily hard-coded for quicker presentation purposes
    discrepancies = ['Variant49196', 'rs188740886', 'chr4:89854340:T:A_ilmnfwd_ilmnF2BT', '6:161575165-TG', 'rs1475032', 'rs6913878']
    # discrepancies = find_unmatched(before_matched)

    snp_name = st.sidebar.selectbox(label = 'SNP Name Choice', label_visibility = 'collapsed', options=discrepancies)

    display_plots(snp_name, before_matched, full_plot = False, before_after = False, title = 'Our Prediction')
    display_plots(snp_name, full_metrics, before_after = False, title = "Illumina's Prediction")

