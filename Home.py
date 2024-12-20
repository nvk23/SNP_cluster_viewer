import numpy as np
import pandas as pd
import streamlit as st
from io import StringIO, BytesIO
from google.cloud import storage

from hold_method import plot_clusters_seaborn, plot_clusters_plotly


st.set_page_config(
    page_title="Array vs. Short-read",
    layout="wide",
    initial_sidebar_state="expanded"
)

# def snp_callback():
#     st.session_state['old_snp_choice'] = st.session_state['snp_choice']
#     st.session_state['snp_choice'] = st.session_state['new_snp_choice']

# Google Cloud file access
# reads in file from google cloud folder
def blob_as_csv(bucket, path, sep='\s+', header='infer'):
    blob = bucket.get_blob(path)
    blob = blob.download_as_bytes()
    blob = str(blob, 'utf-8')
    blob = StringIO(blob)
    df = pd.read_csv(blob, sep=sep, header=header)
    return df

# gets folders from Google Cloud
def get_gcloud_bucket(bucket_name): 
    storage_client = storage.Client(project='gp2-release-terra')
    bucket = storage_client.bucket(bucket_name, user_project='gp2-release-terra')
    return bucket

# pull data from different Google Cloud folders
gp2_data_bucket = get_gcloud_bucket('gp2_working_eu')

st.title('X Chromosome: Array vs. Short-read SNP Comparisons')

par_snps = blob_as_csv(gp2_data_bucket, 'nicole/x_chrom/sampled_100_PAR.csv', sep = ',')
non_par_snps = blob_as_csv(gp2_data_bucket, 'nicole/x_chrom/sampled_100_non_PAR.csv', sep = ',')
master = blob_as_csv(gp2_data_bucket, 'nicole/x_chrom/subset_rel8_IID.csv', sep = ',')
imputed_snps = blob_as_csv(gp2_data_bucket, 'nicole/x_chrom/imputed_variant_key.csv', sep = ',')
wgs_snps = blob_as_csv(gp2_data_bucket, 'nicole/x_chrom/wgs_variant_key.csv', sep = ',')

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
    snp_name = st.selectbox('Non-PAR SNP Selection:', options=non_par_snps.snpID.values)


nba,short_read = st.columns(2)
nba.markdown('### Array Data')
array_options = {'Raw Genotypes': 'raw', 'Imputed Genotypes': 'imputed'}
nba_gt = nba.selectbox('NBA Genotype Selection:', options=['Raw Genotypes', 'Imputed Genotypes'], label_visibility='collapsed')
short_read.markdown('### Short-Read Data')
wgs_gt = short_read.selectbox('Short-read Genotype Selection', options=['Deep Variant Joint Call'], label_visibility='collapsed')

nba_raw_path = f'nicole/x_chrom/nba_{array_options["Raw Genotypes"]}_snps/{snp_name}_main.csv'
nba_raw = blob_as_csv(gp2_data_bucket, nba_raw_path, sep = ',')
nba_raw = nba_raw.merge(master[['IID', 'GP2ID', 'nba_GP2sampleID_r7']], on = 'IID')
gtype_col = 'GT'

wgs_snp_name = wgs_snps[wgs_snps.snpID == snp_name].keep_ID.iloc[0]

app_stop_wgs = False
try:
    # VCF-based path
    # wgs_path = f'nicole/x_chrom/wgs_extracted_snps/{wgs_snp_name}_main.csv'
    # wgs_gtype = 'GT_new'

    # AD-based path
    extra_allele = wgs_snps[wgs_snps.snpID == snp_name].REF.iloc[0]
    wgs_path = f'nicole/x_chrom/wgs_extracted_snps_AD/{wgs_snp_name}_{extra_allele}_main.csv'
    wgs_gtype = 'GT'

    wgs = blob_as_csv(gp2_data_bucket, wgs_path, sep = ',')
    wgs.dropna(inplace = True)
    wgs['snpID'] = snp_name
    wgs = wgs.merge(nba_raw[['snpID', 'GP2ID', 'R', 'Theta']], left_on = ['snpID', 'IID'], right_on = ['snpID', 'GP2ID'], how = 'inner')
except:
    short_read.warning('No WGS data for this SNP')
    app_stop_wgs = True

app_stop_imputed = False
if array_options[nba_gt] == 'imputed':
    imputed_snp_name = imputed_snps[imputed_snps.snpID == snp_name].keep_ID.iloc[0]

    try:
        # VCF-based path
        # imputed_path = f'nicole/x_chrom/imputed_extracted_snps/{imputed_snp_name}_main.csv'
        # gtype_col = 'GT_new'

        # AD-based path
        extra_allele = imputed_snps[imputed_snps.snpID == snp_name].REF.iloc[0]
        imputed_path = f'nicole/x_chrom/imputed_extracted_snps_AD/{imputed_snp_name}_{extra_allele}_main.csv'

        nba_imputed = blob_as_csv(gp2_data_bucket, imputed_path, sep = ',')
        nba_imputed['snpID'] = snp_name

        nba_raw = nba_imputed.merge(nba_raw[['snpID', 'nba_GP2sampleID_r7', 'R', 'Theta']], left_on = ['snpID', 'IID'], right_on = ['snpID', 'nba_GP2sampleID_r7'], how = 'inner')
    except:
        nba.warning('No imputed data for this SNP')
        app_stop_imputed = True

if not app_stop_imputed:
    nba.metric(f'Number of available samples:', "{:.0f}".format(len(np.unique(nba_raw.IID))))
    # nba_plot = plot_clusters_seaborn(nba_raw, x_col='Theta', y_col='R', gtype_col=gtype_col, title = f'{snp_name} NBA {nba_gt}', opacity = 1)['fig']
    # nba.pyplot(nba_plot)
    nba_plot = plot_clusters_plotly(nba_raw, x_col='Theta', y_col='R', gtype_col=gtype_col, title = f'{snp_name} NBA {nba_gt}', opacity = 1)
    nba.plotly_chart(nba_plot)

    # add "compare sex" button at the bottom of plots on both sides
    nba_compare = nba.checkbox('__Compare Biological Sex__', key = 'nba_gender')

    if nba_compare:
        sex1, sex2 = st.columns(2)

        if array_options[nba_gt] == 'imputed':
            male_IID = master[master.biological_sex_for_qc == 'Male'].nba_GP2sampleID_r7.values
        else:
            male_IID = master[master.biological_sex_for_qc == 'Male'].IID.values

        # nba_plot_male = plot_clusters_seaborn(nba_raw[nba_raw.IID.isin(male_IID)], x_col='Theta', y_col='R', gtype_col=gtype_col, title = f'Male {snp_name} NBA {nba_gt}', opacity = 1)['fig']
        # sex1.pyplot(nba_plot_male)
        nba_plot_male = plot_clusters_plotly(nba_raw[nba_raw.IID.isin(male_IID)], x_col='Theta', y_col='R', gtype_col=gtype_col, title = f'Male {snp_name} NBA {nba_gt}', opacity = 1)
        sex1.plotly_chart(nba_plot_male)

        # nba_plot_female = plot_clusters_seaborn(nba_raw[~nba_raw.IID.isin(male_IID)], x_col='Theta', y_col='R', gtype_col=gtype_col, title = f'Female {snp_name} NBA {nba_gt}', opacity = 1)['fig']
        # sex2.pyplot(nba_plot_female)
        nba_plot_female = plot_clusters_plotly(nba_raw[~nba_raw.IID.isin(male_IID)], x_col='Theta', y_col='R', gtype_col=gtype_col, title = f'Female {snp_name} NBA {nba_gt}', opacity = 1)
        sex2.plotly_chart(nba_plot_female)

# Will still plot NBA samples if WGS not available
if app_stop_wgs:
    st.stop()

short_read.metric(f'Number of available samples:', "{:.0f}".format(len(np.unique(wgs.IID))))
# wgs_plot = plot_clusters_seaborn(wgs, x_col='Theta', y_col='R', gtype_col=wgs_gtype, title = f'{snp_name} WGS {wgs_gt}', opacity = 1)['fig']
# short_read.pyplot(wgs_plot)
wgs_plot = plot_clusters_plotly(wgs, x_col='Theta', y_col='R', gtype_col=wgs_gtype, title = f'{snp_name} WGS {wgs_gt}', opacity = 1)
short_read.plotly_chart(wgs_plot)
short_compare = short_read.checkbox('__Compare Biological Sex__', key = 'wgs_gender')

if short_compare:
    sex1, sex2 = st.columns(2)

    male_IID = master[master.biological_sex_for_qc == 'Male'].GP2ID.values

    # wgs_plot_male = plot_clusters_seaborn(wgs[wgs.IID.isin(male_IID)], x_col='Theta', y_col='R', gtype_col=wgs_gtype, title = f'Male {snp_name} WGS {wgs_gt}', opacity = 1)['fig']
    # sex1.pyplot(wgs_plot_male)
    wgs_plot_male = plot_clusters_plotly(wgs[wgs.IID.isin(male_IID)], x_col='Theta', y_col='R', gtype_col=wgs_gtype, title = f'Male {snp_name} WGS {wgs_gt}', opacity = 1)
    sex1.plotly_chart(wgs_plot_male)

    # wgs_plot_female = plot_clusters_seaborn(wgs[~wgs.IID.isin(male_IID)], x_col='Theta', y_col='R', gtype_col=wgs_gtype, title = f'Female {snp_name} WGS {wgs_gt}', opacity = 1)['fig']
    # sex2.pyplot(wgs_plot_female)
    wgs_plot_female = plot_clusters_plotly(wgs[~wgs.IID.isin(male_IID)], x_col='Theta', y_col='R', gtype_col=wgs_gtype, title = f'Female {snp_name} WGS {wgs_gt}', opacity = 1)
    sex2.plotly_chart(wgs_plot_female)