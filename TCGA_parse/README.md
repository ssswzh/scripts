# summarize_compare_hotspot.py

`summarize_compare_hotspot.py` is the upgrade version of `summarize_maf_hotspot.py`. It can summarize the mutation as `summarize_maf_hotspot.py`, and also it can compare two datasets.

Usage:

```
usage: summarize_compare_hotspot.py [-h] --maf MAF --tumor TUMOR --out OUT [--info INFO] [--by {Tumor,Subtype}] [--site {ALL,Primary,Metastasis,Recurrence,Not_Applicable}]
                                    [--variant {ALL,DEL,INS,SNP,DNP,TNP,ONP}]
                                    [--func {ALL,3'Flank,3'UTR,5'Flank,5'UTR,Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Intron,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,RNA,Silent,Splice_Site,Translation_Start_Site}]
                                    [--base BASE] [--mode {Summary,Compare,Both}]

Summarize tumor hotspot according to maf, compare hotspot according to baseline.

Required arguments:
  --maf MAF             maf file
                        if choose to only compare to baseline (--mode==Compare --base), specify summarized file XX.HotSites.tsv, otherwise specify regular maf file
  --tumor TUMOR         tumor type to be summarized, LUAD
  --out OUT             output prefix

Optional arguments:
  --info INFO           sample barcode and disease subtype file
                        if choose --mode==Summary or --mode==Both, specify info file
  --by {Tumor,Subtype}  summarize hotspot by "Tumor" or "Subtype"
  --site {ALL,Primary,Metastasis,Recurrence,Not_Applicable}
                        summarize hotspot by sites in column "SITE", do not specify if you do not have this column
  --variant {ALL,DEL,INS,SNP,DNP,TNP,ONP}
                        Variant_Type type to be reported, separate by comma ',', use 'ALL' for NOT FILTER any type
                        default 'ALL', choose from [ALL,DEL,INS,SNP,DNP,TNP,ONP]
  --func {ALL,3'Flank,3'UTR,5'Flank,5'UTR,Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Intron,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,RNA,Silent,Splice_Site,Translation_Start_Site}
                        Variant_Classification type to be reported, separate by comma ',', use 'ALL' for NOT FILTER any classification
                        default 'Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Intron,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,RNA,Silent,Splice_Site,Translation_Start_Site'
                        choose from [ALL,3'Flank,3'UTR,5'Flank,5'UTR,Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Intron,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,RNA,Silent,Splice_Site,Translation_Start_Site]
  --base BASE           mutation site baseline maf, if specified, compare input maf to baseline maf
                        if --by==subtype, accept key:value pair for baseline maf of each subtype OR one baseline maf for comparing, separate key:value pair by comma ","
                        For types do not have baseline file, use "Other:XXX.tsv" for comparison.
                        if --by==tumor, only accept one baseline maf file
  --mode {Summary,Compare,Both}
                        Summary OR Compare OR Both

Usage:
    Summary:
        summarize_compare_hotspot.py --maf Add.CancerType.maf --info Tumor_subtypes_for_PanCanPathways_9125.txt --tumor BRCA --out BRCA_hotspot/TCGA --by subtype --variant ALL > log
    Compare:
        summarize_compare_hotspot.py --maf TCGA.BRCA.HotSites.tsv --tumor BRCA --out output --by subtype --site Primary --base TCGA.BRCA.HotSites.tsv --mode Compare > log
    Summary and Compare:
        summarize_compare_hotspot.py --maf BRCA/brca_ink4_msk_2021/data_mutation.add_clinical.txt --info BRCA/brca_ink4_msk_2021/data_clinical_sample.info.txt --tumor BRCA --out output --by subtype --site Primary --variant ALL --base TCGA.BRCA.HotSites.tsv --mode Both > log
        summarize_compare_hotspot.py --maf data_mutation.add_clinical.txt --info data_clinical_sample.info.txt --tumor BRCA --out output --by subtype --site Primary --variant ALL --base Her2:BRCA_Her2.HotSites.tsv,LumA:BRCA_LumA.HotSites.tsv,Other:BRCA.HotSites.tsv --mode Both > log

```

The input is MAF format file. `--info` file is tab-delimited file as follow. if choose `--mode==Summary` or `--mode==Both`, info file should be given. If specify `--base`, compare input maf to baseline maf.

|PATIENT_BARCODE    |SAMPLE_BARCODE     |DISEASE    |SUBTYPE        |
|       ---         |        ---        |    ---    |     ---       |
|TCGA-OR-A5J1       |TCGA-OR-A5J1-01    |ACC        |Not_Applicable |
|TCGA-OR-A5J2       |TCGA-OR-A5J2-01    |ACC        |Not_Applicable |
|TCGA-OR-A5J3       |TCGA-OR-A5J3-01    |ACC        |Not_Applicable |
|TCGA-OR-A5J5       |TCGA-OR-A5J5-01    |ACC        |Not_Applicable |


# merge_hotsites.py

`merge_hotsites.py` is used to merge hotspots summarized by `summarize_compare_hotspot.py`. 

The output file can be used as `--base` file in `summarize_compare_hotspot.py`.

Usage:

```
usage: merge_hotsites.py [-h] --files FILES --out OUT [--sort {CoveredSampleSum,AverageCoveredSamplePercent}] [--mode {intersect,union}]

Merge dataframes by columns (like cbind).

Required arguments:
  --files FILES         input files, sepratated by comma ","
  --out OUT             output file name

Optional arguments:
  --sort {CoveredSampleSum,AverageCoveredSamplePercent}
                        sort by condition, default [CoveredSampleSum]
                        CoveredSampleSum: sum of covered sample numbers across datasets
                        AverageCoveredSamplePercent: average of covered sample percent across datasets
  --mode {intersect,union}
                        merge mode, default "union"

Other info:
    Merge by columns:
        ['Chromosome','Start_Position','End_Position','Variant_Classification','Variant_Type','Reference_Allele']
    Rename columns by pattern:
        ['Occurrence','Covered_Sample', 'Covered_Sample_Number', 'Cumulate_Sample_Number', 'Cumulate_Sample_Number_Percent']

Usage:
    python merge_hotsites.py --files TCGA.LUAD.HotSites.tsv,...,data_mutation.Summary.LUAD.HotSites.tsv --out LUAD.merged_baseline.tsv --sort AverageCoveredSamplePercent --mode union > LUAD.merge.log
    
```


compare_hotspot_saturation.r

hotspot_saturation.r

# merge_mutation_clinic.py

```
usage: merge_mutation_clinic.py [-h] --maf MAF --info INFO --out OUT [--col COL]

Merge infomation according to data sample barcode to maf.

Required arguments:
  --maf MAF    maf file
  --info INFO  TCGA baseline, containing several columns
  --out OUT    output file name

Optional arguments:
  --col COL    merge dataframes by column name in info file
```


compare_multiple_hotspot.r



pv_count.py