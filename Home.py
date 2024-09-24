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
from io import StringIO, BytesIO

from hold_method import plot_clusters_seaborn


st.set_page_config(
    page_title="Array vs. Short-read",
    layout="wide",
    initial_sidebar_state="expanded"
)

# def snp_callback():
#     st.session_state['old_snp_choice'] = st.session_state['snp_choice']
#     st.session_state['snp_choice'] = st.session_state['new_snp_choice']


st.title('X Chromosome: Array vs. Short-read SNP Comparisons')

par_snps = pd.read_csv(f'data/sample_info/sampled_100_PAR.csv')
non_par_snps = pd.read_csv(f'data/sample_info/sampled_100_non_PAR.csv')

# Gives user option to choose studies based on their names or target disease
choice = st.radio(
    "Choose if you want to select a PAR or non-PAR SNP",
    ('PAR', 'Non-PAR'), label_visibility='collapsed', horizontal=True) # choose search option
st.write('<style>div.row-widget.stRadio > div{flex-direction:row;justify-content: center;} </style>', unsafe_allow_html=True)

if choice == 'PAR':
    snp_name = st.selectbox('PAR SNP Selection:', options=par_snps.snpID.values)
    st.selectbox(label = 'Non-PAR SNP Selection', disabled = True, options=non_par_snps.snpID.values)
else:
    st.selectbox(label = 'PAR SNP Selection', disabled = True, options=par_snps.snpID.values)
    snp_name = st.selectbox('Non-PAR SNP Selection:', options=par_snps.snpID.values)

nba_raw = pd.read_csv(f'data/nba_raw_snps/{snp_name}_main.csv')


nba,short_read = st.columns(2)
nba.markdown('### Array Data')
nba.selectbox('NBA Genotype Selection:', options=['Raw Genotypes', 'Imputed Genotypes'], label_visibility='collapsed')
# short_read.markdown('### Short-Read Data')

nba.metric(f'Number of available samples:', "{:.0f}".format(len(np.unique(nba_raw.IID))))
# nba.metric(f'Number of available samples:', "{:.0f}".format(len(np.unique(short_read.IID))))

nba_plot = plot_clusters_seaborn(nba_raw, x_col='Theta', y_col='R', gtype_col='GT', title = snp_name, opacity = 1)['fig']
nba.pyplot(nba_plot)