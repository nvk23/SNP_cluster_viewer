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


st.set_page_config(
    page_title="Cluster Buster Evaluation",
    layout="wide",
    initial_sidebar_state="expanded"
)

def plot_clusters(df, x_col='theta', y_col='r', gtype_col='gt', title='snp plot'):
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

def snp_callback():
    st.session_state['old_snp_choice'] = st.session_state['snp_choice']
    st.session_state['snp_choice'] = st.session_state['new_snp_choice']


st.title('Cluster Buster Evaluations')

st.sidebar.markdown('### Choose a MAF Category')
maf_descriptions = st.sidebar.expander("MAF Categories", expanded=False)
with maf_descriptions:
    st.markdown("MAF of 0 is 0 ")
    st.markdown("MAF of 1 is between 0 and 0.005 inclusive")
    st.markdown("MAF of 2 is between 0.005 and 0.01 inclusive")
    st.markdown("MAF of 3 is between 0.01 and 0.05 inclusive")
    st.markdown("MAF of 4 is greater than 0.05")
maf_cat = st.sidebar.selectbox(label = 'Minor Allele Frequency Category Selection', label_visibility = 'collapsed', options=list(range(5)))

checkbox1, checkbox2 = st.columns(2)

if maf_cat == 0:
    prev_nc = checkbox1.checkbox('Show previously NC only')
else:
    prev_nc = False

if not prev_nc:
    metrics = pd.read_csv(f'data/model_052523_maf{maf_cat}')
else:
    metrics = pd.read_csv(f'data/model_052523_maf{maf_cat}_prevNC')
    confid_level = checkbox2.checkbox('Choose by confidence level')

    if confid_level:
        confidence = st.select_slider('Display a confidence level of predictions that is less than:', options=['100', '90', '80', '70', '60'])
        metrics = pd.read_csv(f'data/model_052523_maf{maf_cat}_prevNC_proba{confidence}')

# st.dataframe(metrics)

st.sidebar.markdown('### Choose an individual SNP to display')
snp_name = st.sidebar.selectbox(label = 'Cohort Selection', label_visibility = 'collapsed', options=metrics['snpID'].unique())
full_MAF_plot = st.sidebar.checkbox('Display full MAF plot')
seaborn_plot = st.sidebar.checkbox('Display Seaborn full MAF plot')
full_plots = st.sidebar.checkbox('Display full plots')


metric1,metric2 = st.columns([1,1])
num_snps = len(metrics['snpID'].unique())
num_sample_metrics = len(metrics['GP2sampleID'].unique())

with metric1:
    st.metric(f'Number of available SNPs:', "{:.0f}".format(num_snps))

with metric2:
    st.metric(f'Number of samples with SNP metrics available:',"{:.0f}".format(num_sample_metrics))

snp1 = metrics.loc[metrics['snpID'] == snp_name]

before, after = st.columns(2)

fig_before = plot_clusters(snp1, x_col='Theta', y_col='R', gtype_col='GT', title = 'Before Recluster')['fig']
before.plotly_chart(fig_before, use_container_width=True)

fig_after = plot_clusters(snp1, x_col='Theta', y_col='R', gtype_col='preds_cat', title = 'After Recluster')['fig']
after.plotly_chart(fig_after, use_container_width=True)


title1, title2, title3 = st.columns(3)
maf_full1, maf_full2 = st.columns(2)
tl1, tl2, tl3 = st.columns(3)
full1, full2 = st.columns(2)
color_dict = dict({'NC':'red',
                  'AA':'blue',
                  'BB': 'green',
                  'AB': 'orange'})
if full_MAF_plot:
    title2.markdown('### Full Plots per MAF Category')
    full_metrics = pd.read_csv(f'data/model_052523_maf{maf_cat}_full')
    fig_before = plot_clusters(full_metrics, x_col='Theta', y_col='R', gtype_col='GT', title = 'Before Recluster')['fig']
    maf_full1.plotly_chart(fig_before, use_container_width=True)

    fig_after = plot_clusters(full_metrics, x_col='Theta', y_col='R', gtype_col='preds_cat', title = 'After Recluster')['fig']
    maf_full2.plotly_chart(fig_after, use_container_width=True)

if seaborn_plot:
    fig_maf = plt.figure(figsize=(10, 4))
    sns.scatterplot(x=full_metrics["Theta"], y=full_metrics["R"], hue=full_metrics["GT"], palette = color_dict)
    maf_full1.pyplot(fig_maf, use_container_width=True)

    fig_maf2 = plt.figure(figsize=(10, 4))
    sns.scatterplot(x=full_metrics["Theta"], y=full_metrics["R"], hue=full_metrics["preds_cat"], palette = color_dict)
    maf_full2.pyplot(fig_maf2, use_container_width=True)

if full_plots:
    tl2.markdown('### NC-Only Plot vs. Full Plot')
    nc_metrics = pd.read_csv(f'data/model_052523_nc')
    fig_before = plot_clusters(nc_metrics, x_col='Theta', y_col='R', gtype_col='GT', title = 'Full NC')['fig']
    full1.plotly_chart(fig_before, use_container_width=True)

    full = pd.read_csv(f'data/model_052523_full')
    fig_after = plot_clusters(full, x_col='Theta', y_col='R', gtype_col='GT', title = 'Full Total')['fig']
    full2.plotly_chart(fig_after, use_container_width=True)


#### From GenoTools App - SNP Metrics Page (Additional features if wanted)

# chr_ancestry_select()

# chr_choice = st.session_state['chr_choice']
# ancestry_choice = st.session_state['ancestry_choice']
# selection = f'{ancestry_choice}_{chr_choice}'

# # FIN doesn't have SNP metrics yet (> 50 samples)
# if ancestry_choice ==  'FIN':
#     st.error(f"SNP metrics are not yet available for FIN because less than 50 samples have been released for this ancestry. Please \
#              select a different ancestry!")

# else:
#     metrics_blob_name = f'gp2_snp_metrics/{ancestry_choice}/chr{chr_choice}_metrics.csv'
#     maf_blob_name = f'gp2_snp_metrics/{ancestry_choice}/{ancestry_choice}_maf.afreq'
#     full_maf_blob_name = f'gp2_snp_metrics/full_maf.afreq'

#     if selection not in st.session_state:
#         metrics = pd.read_csv(metrics_blob_name, sep=',')
#         st.session_state[selection] = metrics
#     else:
#         metrics = st.session_state[selection]

#     if f'{ancestry_choice}_maf' not in st.session_state:
#         maf = pd.read_csv(maf_blob_name, sep='\t')
#         st.session_state[f'{ancestry_choice}_maf'] = maf
#     else:
#         maf = st.session_state[f'{ancestry_choice}_maf']

#     if 'full_maf' not in st.session_state:
#         full_maf = pd.read_csv(full_maf_blob_name, sep='\t')
#         st.session_state[f'full_maf'] = full_maf
#     else:
#         full_maf = st.session_state['full_maf']


#     # metrics.columns = ['snpid','r','theta','gentrainscore','gt','chromosome','position','iid','phenotype']

#     metric1,metric2 = st.columns([1,1])

#     num_snps = len(metrics['snpID'].unique())
#     num_sample_metrics = len(metrics['Sample_ID'].unique())

#     with metric1:
#         st.metric(f'Number of available SNPs on Chromosome {chr_choice} for {ancestry_choice}:', "{:.0f}".format(num_snps))

#     with metric2:
#         st.metric(f'Number of {ancestry_choice} samples with SNP metrics available:',"{:.0f}".format(num_sample_metrics))

#     metrics_copy = metrics.copy(deep=True)
#     metrics_copy['snp_label'] = metrics_copy['snpID'] + ' (' + metrics_copy['chromosome'].astype(str) + ':' + metrics_copy['position'].astype(str) + ')'

#     if num_sample_metrics > 0:
#         snp_options = ['Select SNP!']+[snp for snp in metrics_copy['snp_label'].unique()]

#         if 'snp_choice' not in st.session_state:
#             st.session_state['snp_choice'] = snp_options[0]
#         if 'old_snp_choice' not in st.session_state:
#             st.session_state['old_snp_choice'] = ""

#         if st.session_state['snp_choice'] in snp_options:
#             index = snp_options.index(st.session_state['snp_choice'])
        
#         if st.session_state['snp_choice'] not in snp_options:
#             if ((st.session_state['snp_choice'] != 'Select SNP!') and (int(st.session_state['snp_choice'].split('(')[1].split(':')[0]) == st.session_state['old_chr_choice'])):
#                 index = 0
#             else:
#                 st.error(f'SNP: {st.session_state["snp_choice"]} is not availble for {ancestry_choice}. Please choose another SNP!')
#                 index=0

#         st.markdown('### Select SNP for Cluster Plot')

#         st.session_state['snp_choice'] = st.selectbox(label='SNP', label_visibility='collapsed', options=snp_options, index=index, key='new_snp_choice', on_change=snp_callback)

#         if st.session_state['snp_choice'] != 'Select SNP!':
#             snp_df = metrics_copy[metrics_copy['snp_label'] == st.session_state['snp_choice']]
#             snp_df = snp_df.reset_index(drop=True)

#             fig = plot_clusters(snp_df, x_col='Theta', y_col='R', gtype_col='GT', title=st.session_state['snp_choice'])['fig']

#             col1, col2 = st.columns([2.5,1])

#             with col1:
#                 st.plotly_chart(fig, use_container_width=True)

#                 ### FOR EXPORTING CLUSTER PLOT AS PNG
#                 # export_button = st.button('Click here to export displayed SNP metrics cluster plot to .png!')
#                 # if export_button:
#                 #     file_name = f'cluster_plots/{ancestry_choice}_{str(snp_df["snpID"].values[0]).replace(" ","_")}.png'
#                 #     blob = snp_metrics_bucket.blob(file_name)

#                 #     buf = BytesIO()
#                 #     fig.write_image(buf, width=1980, height=1080)
#                 #     blob.upload_from_file(buf, content_type='image/png', rewind=True)

#                 #     st.markdown(f'Cluster plot for {st.session_state["snp_choice"]} written to {snp_metrics_bucket_name}/{file_name}')

#             hide_table_row_index = """<style>thead tr th:first-child {display:none} tbody th {display:none}"""
#             st.markdown(hide_table_row_index, unsafe_allow_html=True)

#             with col2:
#                 st.metric(f'GenTrain Score:', "{:.3f}".format(snp_df['GenTrain_Score'][0]))

#                 within_ancestry_maf = maf[maf['ID'] == snp_df['snpID'].values[0]]
#                 across_ancestry_maf = full_maf[full_maf['ID'] == snp_df['snpID'].values[0]]

#                 st.metric(f'Minor Allele Frequency within {ancestry_choice}:', "{:.3f}".format(within_ancestry_maf['ALT_FREQS'].values[0]))
#                 st.metric(f'Minor Allele Frequency across ancestries:', "{:.3f}".format(across_ancestry_maf['ALT_FREQS'].values[0]))

#                 with st.expander('**Control Genotype Distribution**'):
#                     gt_counts = snp_df[snp_df['phenotype'] == 'Control']['GT'].value_counts().rename_axis('Genotype').reset_index(name='Counts')
#                     gt_rel_counts = snp_df[snp_df['phenotype'] == 'Control']['GT'].value_counts(normalize=True).rename_axis('Genotype').reset_index(name='Frequency')
#                     gt_counts = pd.concat([gt_counts, gt_rel_counts['Frequency']], axis=1)
#                     st.table(gt_counts)

#                 with st.expander('**PD Genotype Distribution**'):
#                     gt_counts = snp_df[snp_df['phenotype'] == 'PD']['GT'].value_counts().rename_axis('Genotype').reset_index(name='Counts')
#                     gt_rel_counts = snp_df[snp_df['phenotype'] == 'PD']['GT'].value_counts(normalize=True).rename_axis('Genotype').reset_index(name='Frequency')
#                     gt_counts = pd.concat([gt_counts, gt_rel_counts['Frequency']], axis=1)
#                     st.table(gt_counts)