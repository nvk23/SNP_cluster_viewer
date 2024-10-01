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
    page_title="SNCA Artifact-Range vs. Not In-Range",
    layout="wide",
    initial_sidebar_state="expanded"
)

# def snp_callback():
#     st.session_state['old_snp_choice'] = st.session_state['snp_choice']
#     st.session_state['snp_choice'] = st.session_state['new_snp_choice']


st.title('SNCA Artifact SNPs vs. Non-Artifact SNPs')

artifact_snps = pd.read_csv(f'data/')
non_artifact_snps = pd.read_csv(f'data/')

artifact,non_artifact = st.columns(2)

# Gives user option to choose studies based on their names or target disease
artifact_name = artifact.selectbox('Artifact SNP Selection:', options=artifact_snps.snpID.values)
non_artifact_name = non_artifact.selectbox('Non-PAR SNP Selection:', options=non_artifact_snps.snpID.values)

artifact.markdown('### Artifact SNPs')
non_artifact.markdown('### Non-Artifact SNPs')

artifact.metric(f'Number of available samples:', "{:.0f}".format(len(np.unique(nba_raw.IID))))
non_artifact.metric(f'Number of available samples:', "{:.0f}".format(len(np.unique(nba_raw.IID))))

art_plot = plot_clusters_seaborn(nba_raw, x_col='Theta', y_col='R', gtype_col='GT', title = f'{snp_name} NBA {nba_gt}', opacity = 1)['fig']
artifact.pyplot(art_plot)