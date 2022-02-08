## dragen-covid-pipeline-epitrax_v2
Getting collection dates and other metadata for EpiTrax upload. 
The script detects the last two files created in IDGenomics_NAS/NextStrain/Epitrax: NGS_Covid_xxxx_xx_xx.csv and All_ncovid_NoNGS_to_xx_xx_xxxx.csv
The main idea: meging tables: covidseq_summary and the NGSxxx, All_ncovidxxx
The scrip needs an argumen: The run name

i.e: Rscript /Volumes/Bioinformatics/. . . ./ . . . ./dragen-covid-pipeline-epitrax_v2.R UT-A01290-220202
the ngs file will be saved at the sumcovidseq folder
