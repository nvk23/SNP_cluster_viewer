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

artifact_snps = pd.read_csv(f'data/sample_info/artifact_SNCA_SNPs.csv')
non_artifact_snps = pd.read_csv(f'data/sample_info/non_artifact_SNCA_SNPs.csv')


st.markdown('*:red[Believed artifact range: 89800000-89870000]*')
artifact,non_artifact = st.columns(2)

artifact.markdown('### Artifact SNPs')
non_artifact.markdown('### Non-Artifact SNPs')

# Gives user option to choose studies based on their names or target disease
artifact_name = artifact.selectbox('Artifact SNP Selection:', options=artifact_snps.snpID.values)
non_artifact_name = non_artifact.selectbox('Non-Artifact SNP Selection:', options=non_artifact_snps.snpID.values)

art_snps = pd.read_csv(f'data/extracted_snps/{artifact_name}_main.csv')
non_art_snps = pd.read_csv(f'data/extracted_snps/{non_artifact_name}_main.csv')

artifact.metric(f'Number of available samples:', len(np.unique(art_snps.IID)))
non_artifact.metric(f'Number of available samples:', len(np.unique(non_art_snps.IID)))

artifact.metric(f'Position on chromosome:', artifact_snps[artifact_snps.snpID == artifact_name].iloc[0].position)
non_artifact.metric(f'Position on chromosome:', non_artifact_snps[non_artifact_snps.snpID == non_artifact_name].iloc[0].position)

art_plot = plot_clusters_seaborn(art_snps, x_col='Theta', y_col='R', gtype_col='GT', title = f'{artifact_name}', opacity = 1)['fig']
artifact.pyplot(art_plot)

non_art_plot = plot_clusters_seaborn(non_art_snps, x_col='Theta', y_col='R', gtype_col='GT', title = f'{non_artifact_name}', opacity = 1)['fig']
non_artifact.pyplot(non_art_plot)

artifact.dataframe(artifact_snps[artifact_snps.snpID == artifact_name])
non_artifact.dataframe(non_artifact_snps[non_artifact_snps.snpID == non_artifact_name])