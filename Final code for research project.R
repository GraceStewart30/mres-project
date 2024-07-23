
installed.packages()
install.packages('BiocManager')
library(BiocManager)
setRepositories()
BiocManager::install('cBioPortalData')  
library(cBioPortalData) 

install.packages('RColorBrewer')
library(RColorBrewer)

install.packages('tidyverse')
library(tidyverse)


#Obtaining data from cBioPortal
cbio = cBioPortal()
studies = getStudies(cbio, buildReport = TRUE)
head(studies)

brca_studies = filter(studies, cancerTypeId =='brca')
brca_tcga = cBioDataPack('brca_tcga_pan_can_atlas_2018')
sampleMap(brca_tcga)

ls('package:cBioPortalData')
experiments(brca_tcga)

brca_tcga_sample_data = colData(brca_tcga)
brca_tcga_data = as_tibble(brca_tcga_sample_data)

zscores = assays(brca_tcga)[['mrna_seq_v2_rsem_zscores_ref_normal_samples']]
View(zscores)
data_summary = lapply(brca_tcga_data, summary)

#heatmaps
BiocManager::install('ComplexHeatmap')
library(ComplexHeatmap)
BiocManager::install('colorRamp2')
library('colorRamp2')

seps = c('DIO1', 'DIO2', 'DIO3', 
         'GPX1', 'GPX2', 'GPX3', 'GPX4', 'GPX6', 
         'SEPX1', 
         'SELENOF', 'C11orf31', 'EPT1', 'SELK', 'SELM', 'SEPN1', 'SELO', 'SELP', 'SELS', 'SELT', 'SELV', 'SEPW1',
         'SEPHS2',
         'TXNRD1', 'TXNRD2', 'TXNRD3IT1')

seps_name_mapping = c('DIO1' = 'DIO1', 'DIO2' = 'DIO2', 'DIO3' = 'DIO3', 
                      'GPX1' = 'GPX1', 'GPX2' = 'GPX2', 'GPX3' = 'GPX3', 'GPX4' = 'GPX4', 'GPX6' = 'GPX6',
                      'SEPX1' = 'MSRB1', 
                      'SELENOF' = 'SELENOF', 'C11orf31' = 'SELENOH', 'EPT1' = 'SELENOI', 'SELK' = 'SELENOK', 'SELM' = 'SELENOM', 'SEPN1' = 'SELENON', 'SELO' = 'SELENOO', 
                      'SELP' = 'SELENOP', 'SELS' = 'SELENOS', 'SELT' = 'SELENOT', 'SELV' = 'SELENOV', 'SEPW1' = 'SELENOW',
                      'SEPHS2' = 'SEPHS2',
                      'TXNRD1' = 'TXNRD1', 'TXNRD2' = 'TXNRD2', 'TXNRD3IT1' = 'TXNRD3')

seps_zscores = zscores[rownames(zscores) %in% seps_for_zscores, ]
rownames(seps_zscores) = seps_name_mapping[rownames(seps_zscores)]
seps_row_distance = dist(seps_zscores, method='euclidean')

subtype_pt_id = data.frame(brca_tcga_data$PATIENT_ID,
                           brca_tcga_data$SUBTYPE)
colnames(subtype_pt_id) = c('PATIENT_ID', 'SUBTYPE')
subtype_anno = subtype_pt_id
rownames(subtype_anno) = subtype_anno[,1]
rownames(subtype_anno) = paste0(rownames(subtype_anno), '-01')   
subtype_anno[,1] = NULL
subtype_anno = t(subtype_anno)
subtype_anno = subtype_anno[,colnames(subtype_anno) %in% colnames(seps_zscores)] 
subtype_anno = as.data.frame(subtype_anno)
colnames(subtype_anno) = c('SUBTYPE')
subtype_anno$SUBTYPE = as.factor(subtype_anno$SUBTYPE)

#subtype and tumour type
column_anno_age_tt = column_anno
column_anno_age_tt$TUMOUR_TYPE = brca_tcga_data$TUMOR_TYPE[brca_tcga_data$SAMPLE_ID %in% row.names(column_anno)]
column_anno_age_tt$TUMOUR_TYPE = as.factor(column_anno_age_tt$TUMOUR_TYPE)

subtypes = levels(subtype_anno$SUBTYPE)
subtypes = as.data.frame(subtypes)
subtype_palette = brewer.pal(5, 'Dark2')
tumour_type = levels(column_anno_age_tt$TUMOUR_TYPE)
tumour_type = as.data.frame(tumour_type)
tumour_type_palette = brewer.pal(9, 'Paired')

age_tt_sum = lapply(column_anno_age_tt, summary)

age_tt_anno = HeatmapAnnotation('TUMOUR TYPE' = column_anno_age_tt$TUMOUR_TYPE,
                                'AGE' = column_anno_age_tt$AGE,
                                'SUBTYPE' = column_anno_age_tt$SUBTYPE,
                                col = list('SUBTYPE' = setNames(subtype_palette, subtypes$subtypes),
                                           'AGE' = colorRamp2(seq(min(column_anno_age_tt$AGE), max(column_anno_age_tt$AGE), length.out = length(age_pal)), age_pal),
                                           'TUMOUR TYPE' = setNames(tumour_type_palette, tumour_type$tumour_type)))

st_tt_anno = age_tt_anno = HeatmapAnnotation('TUMOUR TYPE' = column_anno_age_tt$TUMOUR_TYPE,
                                             'SUBTYPE' = column_anno_age_tt$SUBTYPE,
                                             col = list('SUBTYPE' = setNames(subtype_palette, subtypes$subtypes),
                                                        'TUMOUR TYPE' = setNames(tumour_type_palette, tumour_type$tumour_type)))


#Hypoxia
column_anno_hypoxia = column_anno
column_anno_hypoxia$BUFFA_HYPOXIA_SCORE = brca_tcga_data$BUFFA_HYPOXIA_SCORE[brca_tcga_data$SAMPLE_ID %in% row.names(column_anno)]
column_anno_hypoxia$RAGNUM_HYPOXIA_SCORE = brca_tcga_data$RAGNUM_HYPOXIA_SCORE[brca_tcga_data$SAMPLE_ID %in% row.names(column_anno)]

hypoxia_sum = lapply(column_anno_hypoxia, summary)

buffa_pal = brewer.pal(11, 'PiYG')
ragnum_pal = brewer.pal(11, 'RdBu')

hypoxia_anno = HeatmapAnnotation('RAGNUM SCORE' = column_anno_hypoxia$RAGNUM_HYPOXIA_SCORE,
                                 'BUFFA SCORE' = column_anno_hypoxia$BUFFA_HYPOXIA_SCORE,
                                 'SUBTYPE' = column_anno_hypoxia$SUBTYPE,
                                 col = list('SUBTYPE' = setNames(subtype_palette, subtypes$subtypes),
                                            'BUFFA SCORE' = colorRamp2(seq(min(column_anno_hypoxia$BUFFA_HYPOXIA_SCORE, na.rm = TRUE), max(column_anno_hypoxia$BUFFA_HYPOXIA_SCORE, na.rm = TRUE), length.out = length(buffa_pal)), buffa_pal),
                                            'RAGNUM SCORE' = colorRamp2(seq(min(column_anno_hypoxia$RAGNUM_HYPOXIA_SCORE, na.rm = TRUE), max(column_anno_hypoxia$RAGNUM_HYPOXIA_SCORE, na.rm = TRUE), length.out = length(ragnum_pal)), ragnum_pal)))


#Seps heatmap
heatmap_cols = colorRamp2(seq(from=-5, to=5, length.out=101), colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(101))

Heatmap(seps_zscores,
        col = heatmap_cols,
        name = 'Z SCORE',
        cluster_rows = TRUE,
        clustering_distance_rows = seps_row_distance,
        top_annotation = st_tt_anno,
        cluster_columns = TRUE,
        show_column_names = FALSE)

#sep synthesis heatmap

sep_synth = c('SECISBP2',
              'EEFSEC',
              'RPL30',
              'RPL7',
              'SCLY',
              'SEPHS1',
              'SEPHS2',
              'PSTK',
              'SEPSECS',
              'SELP')
sep_synth_namemapping = c('SECISBP2' ='SECISBP2',
                          'EEFSEC' = 'EEFSEC',
                          'RPL30' = 'RPL30',
                          'RPL7' = 'RPL7',
                          'SCLY' = 'SCLY',
                          'SEPHS1' = 'SEPHS1',
                          'SEPHS2' = 'SEPHS2',
                          'PSTK' = 'PSTK', 
                          'SEPSECS' = 'SEPSECS',
                          'SELP' = 'SELENOP')

sep_synth_zscores = zscores[rownames(zscores) %in% sep_synth, ]
rownames(sep_synth_zscores) = sep_synth_namemapping[rownames(sep_synth_zscores)]
sep_synth_row_distance = dist(sep_synth_zscores, method='euclidean')

Heatmap(sep_synth_zscores,
        col = heatmap_cols,
        name = 'Z SCORE',
        cluster_rows = TRUE,
        clustering_distance_rows = sep_synth_row_distance,
        cluster_columns = TRUE,
        top_annotation = st_tt_anno,
        show_column_names = FALSE)

#seps and synth heatmap

seps_and_synth = c('DIO1', 'DIO2', 'DIO3', 
                   'GPX1', 'GPX2', 'GPX3', 'GPX4', 'GPX6', 
                   'SEPX1', 
                   'SELENOF', 'C11orf31', 'EPT1', 'SELK', 'SELM', 'SEPN1', 'SELO', 'SELP', 'SELS', 'SELT', 'SELV', 'SEPW1',
                   'SEPHS2',
                   'TXNRD1', 'TXNRD2', 'TXNRD3IT1',
                   'SECISBP2', 'EEFSEC', 'RPL30', 'RPL7', 'SCLY', 'SEPHS1', 'PSTK', 'SEPSECS')
seps_and_synth_namemapping = c('DIO1' = 'DIO1', 'DIO2' = 'DIO2', 'DIO3' = 'DIO3', 
                               'GPX1' = 'GPX1', 'GPX2' = 'GPX2', 'GPX3' = 'GPX3', 'GPX4' = 'GPX4', 'GPX6' = 'GPX6',
                               'SEPX1' = 'MSRB1', 
                               'SELENOF' = 'SELENOF', 'C11orf31' = 'SELENOH', 'EPT1' = 'SELENOI', 'SELK' = 'SELENOK', 'SELM' = 'SELENOM', 'SEPN1' = 'SELENON', 'SELO' = 'SELENOO', 
                               'SELP' = 'SELENOP', 'SELS' = 'SELENOS', 'SELT' = 'SELENOT', 'SELV' = 'SELENOV', 'SEPW1' = 'SELENOW',
                               'SEPHS2' = 'SEPHS2',
                               'TXNRD1' = 'TXNRD1', 'TXNRD2' = 'TXNRD2', 'TXNRD3IT1' = 'TXNRD3',
                               'SECISBP2' = 'SECISBP2', 'EEFSEC' = 'EEFSEC', 'RPL30' = 'RPL30', 'RPL7' = 'RPL7', 
                               'SCLY' = 'SCLY', 'SEPHS1' = 'SEPHS1', 'PSTK' = 'PSTK', 'SEPSECS' = 'SEPSECS')

seps_and_synth_zscores = zscores[rownames(zscores) %in% seps_and_synth, ]
rownames(seps_and_synth_zscores) = seps_and_synth_namemapping[rownames(seps_and_synth_zscores)]
seps_and_synth_row_distance = dist(seps_and_synth_zscores, method = 'euclidean')


Heatmap(seps_and_synth_zscores,
        col = heatmap_cols,
        name = 'Z SCORE',
        cluster_rows = TRUE,
        clustering_distance_rows = seps_and_synth_row_distance,
        top_annotation = st_tt_anno,
        show_column_names = FALSE)

#methylation

methylation = metadata(brca_tcga)[['methylation_hm27_hm450_merged']]
methylation_table = as.tibble(methylation)
rownames(methylation_data) = methylation_data$NAME

seps_and_synth_for_methylation = c('DIO1', 'DIO2', 'DIO3', 
                                   'GPX1', 'GPX2', 'GPX3', 'GPX4', 
                                   'SEPX1', 
                                   'HS2ST1;SEP15', 'SEPN1', 'SELP', 'SEPW1',
                                   'SEPHS2',
                                   'TXNRD1', 'TXNRD2;COMT',
                                   'SECISBP2', 'RPL7', 'RPL7;RDH10', 'SEPHS1', 'PSTK', 'SEPSECS')

seps_and_synth_methylation_namemapping = c('DIO1' = 'DIO1', 'DIO2' = 'DIO2', 'DIO3' = 'DIO3', 
                                           'GPX1' = 'GPX1', 'GPX2' = 'GPX2', 'GPX3' = 'GPX3', 'GPX4' = 'GPX4', 
                                           'SEPX1' = 'MSRB1', 
                                           'HS2ST1;SEP15' = 'HS2ST1;SELENOF', 'SEPN1' = 'SELENON', 'SELP' = 'SELENOP', 'SEPW1' = 'SELENOW',
                                           'SEPHS2' = 'SEPHS2',
                                           'TXNRD1' = 'TXNRD1', 'TXNRD2;COMT' = 'TXNRD2;COMT', 
                                           'SECISBP2' = 'SECISBP2', 'RPL7' = 'RPL7', 'RPL7;RDH10' = 'RPL7;RDH10',
                                           'SEPHS1' = 'SEPHS1', 'PSTK' = 'PSTK', 'SEPSECS' = 'SEPSECS')

seps_and_synth_methylation = methylation_data[methylation_data$NAME %in% seps_and_synth_for_methylation,]
seps_and_synth_methylation_matrix = seps_and_synth_methylation[, -c(2, 3, 4)]
seps_and_synth_methylation_matrix = as.matrix(seps_and_synth_methylation_matrix)
rownames(seps_and_synth_methylation_matrix) = seps_and_synth_methylation$ENTITY_STABLE_ID
seps_and_synth_methylation_matrix = seps_and_synth_methylation_matrix[, -c(1)]
na_positions_methylation = which(is.na(seps_and_synth_methylation_matrix), arr.ind = TRUE)
seps_and_synth_methylation_distance = dist(seps_and_synth_methylation_matrix, method = 'euclidean')
mode(seps_and_synth_methylation_matrix) = 'numeric'

methylation_cols = colorRamp2(seq(from=0, to=1, length.out=101), colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(101))

entity_ID_name_sep_and_synth = seps_and_synth_methylation[, c(1, 2)]
entity_ID_name_sep_and_synth = as.data.frame(entity_ID_name_sep_and_synth)

entity_ID_name_sep_and_synth$GENE_NAME = seps_and_synth_methylation_namemapping[entity_ID_name_sep_and_synth$NAME]
entity_ID_name_sep_and_synth$GENE_NAME = as.factor(entity_ID_name_sep_and_synth$GENE_NAME)
entity_ID_name_sep_and_synth$CpG_LOCATION = seps_and_synth_methylation$DESCRIPTION
entity_ID_name_sep_and_synth$CpG_LOCATION = as.factor(entity_ID_name_sep_and_synth$CpG_LOCATION)

library(pals)
cpg_pal = kelly(13)

names(cpg_pal) = unique(entity_ID_name_sep_and_synth$CpG_LOCATION)
cpg_pal[entity_ID_name_sep_and_synth$CpG_LOCATION]

left_anno = rowAnnotation('CpG LOCATION' = entity_ID_name_sep_and_synth$CpG_LOCATION,
                          show_annotation_name = FALSE,
                          col = list('CpG LOCATION' = cpg_pal[entity_ID_name_sep_and_synth$CpG_LOCATION]))

Heatmap(seps_and_synth_methylation_matrix,
        col = methylation_cols,
        name = 'METHYLATION',
        show_column_names = FALSE,
        show_row_names = FALSE,
        cluster_rows = F,
        right_annotation = left_anno)
seps_and_synth_for_zscores_and_methylation = c('DIO1', 'DIO2', 'DIO3', 
                                               'GPX1', 'GPX2', 'GPX3', 'GPX4', 
                                               'MSRB1', 
                                               'SELENOF', 'SELENON', 'SELENOP', 'SELENOW',
                                               'SEPHS2',
                                               'TXNRD1', 'TXNRD2',
                                               'SECISBP2', 'RPL7', 'SEPHS1', 'PSTK', 'SEPSECS')

methylation_only_seps_zscores = zscores_for_methylation[rownames(zscores_for_methylation) %in% seps_and_synth_for_zscores_and_methylation,]
methylation_only_seps_zscores_ordered = methylation_only_seps_zscores[order(rownames(methylation_only_seps_zscores)),]

seps_methylation = methylation_data[methylation_data$NAME %in% seps_for_methylation,]
seps_methylation_matrix = seps_methylation[, -c(2, 3, 4)]
seps_methylation_matrix = as.matrix(seps_methylation_matrix)
rownames(seps_methylation_matrix) = seps_methylation$ENTITY_STABLE_ID
seps_methylation_matrix = seps_methylation_matrix[, -c(1)]
na_positions_methylation = which(is.na(seps_methylation_matrix), arr.ind = TRUE)
seps_methylation_distance = dist(seps_methylation_matrix, method = 'euclidean')
mode(seps_methylation_matrix) = 'numeric'


entity_ID_name_seps = seps_methylation[, c(1, 2)]
entity_ID_name_seps = as.data.frame(entity_ID_name_seps)
rownames(entity_ID_name_seps) = entity_ID_name_seps$ENTITY_STABLE_ID
as.factor(entity_ID_name_seps$NAME)

seps_methylation_row_anno = rowAnnotation('GENE NAME' = entity_ID_name_seps$NAME,
                                          col = list('GENE NAME' = setNames(seps_gene_name_pal, entity_ID_name_seps$NAME)))

samples_column_names = colnames(seps_and_synth_methylation_matrix)
samples_column_names = gsub('\\.', '-', samples_column_names)
missing_columns = setdiff(samples_column_names, colnames(seps_and_synth_zscores))
samples_column_names = setdiff(samples_column_names, missing_columns)
zscores_for_methylation = seps_and_synth_zscores[, samples_column_names, drop = FALSE]

colnames(seps_and_synth_methylation_matrix) = gsub('\\.', '-', colnames(seps_and_synth_methylation_matrix))
methylation_for_zscores = seps_and_synth_methylation_matrix[ ,colnames(seps_and_synth_methylation_matrix) %in% samples_column_names]


zscores_for_methylation

identical(colnames(zscores_for_methylation), colnames(methylation_for_zscores))

#zscores heatmap

#time to make some annos :)
zscores_and_methylation_anno = HeatmapAnnotation('TUMOUR TYPE' = methylation_anno_st_tt$TUMOUR_TYPE,
                                                 'SUBTYPE' =  methylation_anno_st_tt$SUBTYPE,
                                                 col = list('SUBTYPE' = setNames(subtype_palette, subtypes$subtypes),
                                                            'TUMOUR TYPE' = setNames(tumour_type_palette, tumour_type$tumour_type)))


zscores_and_methylation_anno_ordered = HeatmapAnnotation('TUMOUR TYPE' = methylation_anno_st_ordered$TUMOUR_TYPE,
                                                         'SUBTYPE' =  methylation_anno_st_ordered$SUBTYPE,
                                                         col = list('SUBTYPE' = setNames(subtype_palette, subtypes$subtypes),
                                                                    'TUMOUR TYPE' = setNames(tumour_type_palette, tumour_type$tumour_type)))



zscores_heatmap = 
  Heatmap(methylation_only_seps_zscores_ordered,
          name = 'Z SCORES',
          col = heatmap_cols,
          top_annotation = zscores_and_methylation_anno,
          cluster_rows = FALSE,
          show_column_names = FALSE)

methylation_heatmap =
Heatmap(methylation_for_zscores,
        name = 'METHYLATION',
        col = methylation_cols,
        left_annotation = left_anno,
        cluster_rows = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE)


zscores_heatmap %v% methylation_heatmap
