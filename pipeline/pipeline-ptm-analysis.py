#################################################################
#################################################################
###############  ################
#################################################################
#################################################################
##### Author: Denis Torre
##### Affiliation: Ma'ayan Laboratory,
##### Icahn School of Medicine at Mount Sinai

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
from ruffus import *
import sys, glob
import pandas as pd
import numpy as np
import rpy2.robjects as robjects

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
sys.path.append('pipeline/scripts')
import Support as S
import PipelinePtmAnalysis as P

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####
rawData = 'rawdata/lim_log2_ratio_data.txt'

##### 2. R Connection #####
rSource = 'pipeline/scripts/pipeline-ptm-analysis.R'
r = robjects.r
r.source(rSource)

#######################################################
#######################################################
########## S1. Process Data
#######################################################
#######################################################

#############################################
########## 1. Process Data
#############################################

@follows(mkdir('f1-data.dir'))

@subdivide(rawData,
	   	   formatter(),
	   	   'f1-data.dir/-*ptm.txt',
	   	   'f1-data.dir')

def processData(infile, outfiles, outfileRoot):

	# Get raw dataframe
	raw_dataframe = pd.read_table(infile)

	# Name index
	raw_dataframe.index.name = 'ptm'

	# Get untreated columns
	untreated_columns = [x for x in raw_dataframe.columns if 'Channel' in x]

	# Get untreated dataframe
	dataframes = {}
	dataframes['untreated'] = raw_dataframe[untreated_columns]

	# Get treated dataframe
	dataframes['treated'] = raw_dataframe[list(set(raw_dataframe.columns)-set(untreated_columns))]

	# Write
	for key, dataframe in dataframes.iteritems():
		dataframe.to_csv('{outfileRoot}/{key}-ptm.txt'.format(**locals()), index=True, sep='\t')

#######################################################
#######################################################
########## S2. Sample Metadata
#######################################################
#######################################################

#############################################
########## 1. Get Treated Metadata
#############################################

@follows(mkdir('f2-metadata.dir'))

@files('f1-data.dir/treated-ptm.txt',
	   'f2-metadata.dir/treated-sample_metadata.txt')

def getTreatedMetadata(infile, outfile):

	# Read dataframe
	treated_dataframe = pd.read_table(infile, index_col='ptm')

	# Get metadata dataframe
	metadata_dataframe = pd.DataFrame({x: {'cell_line': x.split('.')[0], 'timepoint': x.split('.')[1] if len(x.split('.'))>3 else None, 'concentration': x.split('.')[2] if len(x.split('.'))>3 else None, 'drug': x.split('.')[3] if len(x.split('.'))>3 else x.split('.')[1]} for x in treated_dataframe.columns}).T
	metadata_dataframe.index.name = 'sample'

	# Replace concentration
	metadata_dataframe['concentration'] = [x.replace('1uM', '1').replace('100nM', '0.1') if isinstance(x, str) else x for x in metadata_dataframe['concentration']]
	metadata_dataframe['timepoint'] = [x.replace('hr', '') if isinstance(x, str) else x for x in metadata_dataframe['timepoint']]
	metadata_dataframe['drug'] = [x.replace('Control', 'control').replace('control1', 'control').replace('control2', 'control').replace('gleevec1', 'gleevec').replace('gleevec2', 'gleevec').replace('gleevec3', 'gleevec') if isinstance(x, str) else x for x in metadata_dataframe['drug']]
	
	# Write
	metadata_dataframe.to_csv(outfile, sep='\t')

#############################################
########## 2. Get Untreated Metadata
#############################################

@files('f1-data.dir/untreated-ptm.txt',
	   'f2-metadata.dir/untreated-sample_metadata.txt')

def getUnreatedMetadata(infile, outfile):

	# Read dataframe
	untreated_dataframe = pd.read_table(infile, index_col='ptm')

	# Get metadata dataframe
	metadata_dataframe = pd.DataFrame({x: {'cell_line': x.split('.', 3)[-1], 'channel': x.split('.')[1]} for x in untreated_dataframe.columns}).T
	metadata_dataframe.index.name = 'sample'
	
	# Write
	metadata_dataframe.to_csv(outfile, sep='\t')

#############################################
########## 2. Get Modification Type
#############################################

@files('f1-data.dir/treated-ptm.txt',
	   'f2-metadata.dir/ptm-metadata.txt')

def getPtmMetadata(infile, outfile):

	# Read dataframe
	treated_dataframe = pd.read_table(infile, index_col='ptm')

	# Modification dataframe
	modification_dataframe = pd.DataFrame({x: {'modification': x.split(' ')[1], 'gene_symbol': x.split(' ')[0]} for x in treated_dataframe.index}).T

	# Add type
	modification_type_dict = {'ack': 'acetylation', 'p': 'phosphorylation', 'kme': 'methylation', 'rme': 'methylation'}
	modification_dataframe['modification_type'] = [modification_type_dict[x] for x in modification_dataframe['modification']]
	modification_dataframe.index.name = 'ptm'

	# Write
	modification_dataframe.to_csv(outfile, sep='\t')

#############################################
########## 3. Get PTM Counts
#############################################

@transform(getPtmMetadata,
		   suffix('metadata.txt'),
		   'gene_counts.txt')

def getPtmCounts(infile, outfile):

	# Read dataframe
	modification_dataframe = pd.read_table(infile, index_col='ptm')

	# Get counts
	ptm_count_dataframe = modification_dataframe.groupby(['gene_symbol', 'modification_type']).size().rename('count').to_frame().reset_index()

	# Save counts
	ptm_count_dataframe.to_csv(outfile, sep='\t')

#######################################################
#######################################################
########## S3. PTM Correlation
#######################################################
#######################################################

#############################################
########## 1. Get Correlations
#############################################

@follows(mkdir('f3-ptm_correlations.dir'))

@transform(glob.glob('f1-data.dir/*'),
		   regex(r'.*/(.*)-ptm.txt'),
		   r'f3-ptm_correlations.dir/\1-ptm_correlations_top.txt')

def getPtmCorrelations(infile, outfile):

	# Read dataframe
	ptm_dataframe = pd.read_table(infile, index_col='ptm')

	# Get top genes
	topGenes = ptm_dataframe.apply(np.var, axis=1).sort_values(ascending=False).index[:5000]

	# Get correlations
	print 'getting correlations...'
	correlation_dataframe = ptm_dataframe.loc[topGenes].T.corr(method='spearman')
	np.fill_diagonal(correlation_dataframe.values, np.nan)

	# Melt
	print 'melting...'
	correlation_dataframe_melted = pd.melt(correlation_dataframe.reset_index(), id_vars='ptm', var_name='target_ptm', value_name='corr').dropna()

	# Save
	correlation_dataframe_melted.to_csv(outfile, sep='\t', index=False)

#############################################
########## 2. Annotated Correlations
#############################################

@transform(glob.glob('f3-ptm_correlations.dir/*'),
		   suffix('.txt'),
		   add_inputs(getPtmMetadata),
		   '_annotated.txt')

def annotatePtmCorrelations(infiles, outfile):

	# Split infiles
	correlation_infile, annotation_infile = infiles

	# Read dataframe
	correlation_dataframe_melted = pd.read_table(correlation_infile)
	ptm_annotation_dataframe = pd.read_table(annotation_infile, index_col='ptm')

	# Annotate
	annotated_dataframe = correlation_dataframe_melted.merge(ptm_annotation_dataframe, left_on='ptm', right_index=True).rename(columns={'gene_symbol': 'source_gene_symbol', 'modification_type': 'source_modification_type'}).merge(ptm_annotation_dataframe, left_on='target_ptm', right_index=True).rename(columns={'gene_symbol': 'target_gene_symbol', 'modification_type': 'target_modification_type'})[['ptm', 'target_ptm', 'source_gene_symbol', 'target_gene_symbol', 'source_modification_type', 'target_modification_type', 'corr']].rename(columns={'ptm': 'source_ptm'})
	annotated_dataframe['same_gene'] = [x==y for x, y in annotated_dataframe[['source_gene_symbol', 'target_gene_symbol']].as_matrix()]
	
	# Write
	annotated_dataframe.to_csv(outfile, sep='\t', index=False)

##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=2, verbose=1)
print('Done!')
