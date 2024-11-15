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
master = pd.read_csv(f'data/sample_info/subset_rel8_IID.csv')

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


nba,short_read = st.columns(2)
nba.markdown('### Array Data')
array_options = {'Raw Genotypes': 'raw', 'Imputed Genotypes': 'imputed'}
nba_gt = nba.selectbox('NBA Genotype Selection:', options=['Raw Genotypes', 'Imputed Genotypes'], label_visibility='collapsed')
# short_read.markdown('### Short-Read Data')

nba_raw = pd.read_csv(f'data/nba_{array_options[nba_gt]}_snps/{snp_name}_main.csv')
nba.metric(f'Number of available samples:', "{:.0f}".format(len(np.unique(nba_raw.IID))))
# nba.metric(f'Number of available samples:', "{:.0f}".format(len(np.unique(short_read.IID))))

nba_plot = plot_clusters_seaborn(nba_raw, x_col='Theta', y_col='R', gtype_col='GT', title = f'{snp_name} NBA {nba_gt}', opacity = 1)['fig']
nba.pyplot(nba_plot)

# add "compare sex" button at the bottom of plots on both sides
nba_compare = nba.checkbox('Compare Biological Sex')
# short_compare = short_read.radio('Compare Biological Sex')

if nba_compare:
    sex1, sex2 = st.columns(2)
    male_IID = master[master.biological_sex_for_qc == 'Male'].IID.values
    nba_plot_male = plot_clusters_seaborn(nba_raw[nba_raw.IID.isin(male_IID)], x_col='Theta', y_col='R', gtype_col='GT', title = f'Male {snp_name} NBA {nba_gt}', opacity = 1)['fig']
    sex1.pyplot(nba_plot_male)

    nba_plot_female = plot_clusters_seaborn(nba_raw[~nba_raw.IID.isin(male_IID)], x_col='Theta', y_col='R', gtype_col='GT', title = f'Female {snp_name} NBA {nba_gt}', opacity = 1)['fig']
    sex2.pyplot(nba_plot_female)