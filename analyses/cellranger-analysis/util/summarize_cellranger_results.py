#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author:
	Name: 			Antonia Chroni
	Email: 			antonia.chroni@stjude.org
	Affiliation: 	St. Jude Children's Research Hospital, Memphis, TN
	Date: 			April 10th, 2025
"""


import os, sys, argparse, glob, numpy, pandas


def dir_path(string):
	if os.path.isdir(string):
		return string
	else:
		raise NotADirectoryError(string)

#Beginning the argparse module to extract command-line arguments
parser = argparse.ArgumentParser(description="This is a script that will summarize cellranger count results from at least one cellranger output and create a summary within the 4_reports directory of the project. It accepts a data directory that contains at least one cellranger count results.")
#Creating the following optional arguments; some have default values
parser.add_argument('--dir', type=dir_path, help='Data directory path that contains individually named cellranger count results for samples', required=True)
parser.add_argument('--outdir', type=dir_path, help='Create all output files in the specified output directory. Please note that this directory must exist as the program will not create it.', required=True)
parser.add_argument('--genome', type=str, help='Only specify the genome you want to recover data from a multiple genome alignment', required=False)

#Converts argument strings to objects and assigns them as attributes of the namespace; e.g. --id -> args.id
args = parser.parse_args()


MasterDF = pandas.DataFrame()

for filename in glob.glob(os.path.join(args.dir, "*", "outs", "summary.csv")):
    # print(filename)
    
    df = pandas.read_csv(filename)
    df = df.replace(",", "", regex=True)
    df = df.replace("%", "", regex=True)
    df = df.astype('float')
    
    SampleID = filename.split("/outs")[0].split("/")[-1]
#     print(SampleID)
    df["Sample ID"] = SampleID
    
    
    Warnings = ""
    MajorWarnings = ""
    TotalWarnings = 0
    
    
    # Sequencing Metrics 
    if df.iloc[0]["Valid Barcodes"] < 75:
        Warnings = Warnings + "Valid Barcodes < 75%, "
        TotalWarnings += 1
        
    if df.iloc[0]["Q30 Bases in Barcode"] < 65:
        Warnings = Warnings + "Q30 Bases in Barcode < 65%, "
        TotalWarnings += 1
    
    if df.iloc[0]["Q30 bases in read 1"] < 65:
        Warnings = Warnings + "Q30 Bases in RNA Read < 65%, "
        TotalWarnings += 1
        
    if df.iloc[0]["Q30 bases in read 2"] < 65:
        Warnings = Warnings + "Q30 Bases in UMI < 65%, "
        TotalWarnings += 1
    
    if df.iloc[0]["Q30 bases in sample index i1"] < 90:
        Warnings = Warnings + "Q30 Bases in UMI < 90%, "
        TotalWarnings += 1

# Cell Metrics
if args.genome == None:
        if 500 < df.iloc[0]["Estimated Number of Cells"] < 10000:
            Warnings = Warnings + "Estimated Number of Cells < 500 or > 10000, "
            TotalWarnings += 1
elif args.genome == "GRCh38":
        if 500 < df.iloc[0]["GRCh38 Estimated Number of Cells"] < 10000:
            Warnings = Warnings + "GRCh38 Estimated Number of Cells < 500 or > 10000, "
            TotalWarnings += 1
elif args.genome == "GRCm39":
        if 500 < df.iloc[0]["GRCm39 Estimated Number of Cells"] < 10000:
            Warnings = Warnings + "GRCm39 Estimated Number of Cells < 500 or > 10000, "
            TotalWarnings += 1

if args.genome == None:
        if df.iloc[0]["Estimated Number of Cells"] < 100:
            MajorWarnings = MajorWarnings + "Estimated Number of Cells < 100, "
            TotalWarnings += 1
elif args.genome == "GRCh38":
        if df.iloc[0]["GRCh38 Estimated Number of Cell Partitions"] < 100:
            MajorWarnings = MajorWarnings + "GRCh38 Estimated Number of Cell Partitions < 100, "
            TotalWarnings += 1
elif args.genome == "GRCm39":
        if df.iloc[0]["GRCm39 Estimated Number of Cell Partitions"] < 100:
            MajorWarnings = MajorWarnings + "GRCm39 Estimated Number of Cell Partitions < 100, "
            TotalWarnings += 1

if args.genome == None:
        if df.iloc[0]["Fraction of high-quality fragments in cells"] < 40:
            MajorWarnings = MajorWarnings + "Fraction of high-quality fragments in cells < 40%, "
            TotalWarnings += 1
elif args.genome == "GRCh38":
        if df.iloc[0]["GRCh38 Fraction of high-quality fragments in cells"] < 40:
            MajorWarnings = MajorWarnings + "GRCh38 Fraction of high-quality fragments in cells < 40%, "
            TotalWarnings += 1
elif args.genome == "GRCm39":
        if df.iloc[0]["GRCm39 Fraction of high-quality fragments in cells"] < 40:
            MajorWarnings = MajorWarnings + "GRCm39 Fraction of high-quality fragments in cells < 40%, "
            TotalWarnings += 1
    
if args.genome == None:
        if df.iloc[0]["Fraction of transposition events in peaks in cells"] < 15:
            Warnings = Warnings + "Fraction of transposition events in peaks in cells < 15%, "
            TotalWarnings += 1
elif args.genome == "GRCh38":
        if df.iloc[0]["GRCh38 Fraction of transposition events in peaks in cells"] < 15:
            Warnings = Warnings + "GRCh38 Fraction of transposition events in peaks in cells < 15%, "
            TotalWarnings += 1
elif args.genome == "GRCm39":
        if df.iloc[0]["GRCm39 Fraction of transposition events in peaks in cells"] < 15:
            Warnings = Warnings + "GRCm39 Fraction of transposition events in peaks in cells < 15%, "
            TotalWarnings += 1
    
    
# Library Complexity Metric
if args.genome == None:
        if df.iloc[0]["Percent duplicates"] < 30:
            MajorWarnings = MajorWarnings + "Percent duplicates < 30%, "
            TotalWarnings += 1
elif args.genome == "GRCh38":
        if df.iloc[0]["GRCh38 Percent duplicates"] < 30:
            MajorWarnings = MajorWarnings + "GRCh38 Percent duplicates < 30%, "
            TotalWarnings += 1
elif args.genome == "GRCm39":
        if df.iloc[0]["GRCm39 Percent duplicates"] < 30:
            MajorWarnings = MajorWarnings + "GRCm39 Percent duplicates < 30%, "
            TotalWarnings += 1
    
    
# Mapping Metrics
if args.genome == None:
        if df.iloc[0]["Confidently mapped read pairs"] < 80:
            Warnings = Warnings + "Confidently mapped read pairs < 80%, "
            TotalWarnings += 1
elif args.genome == "GRCh38":
        if df.iloc[0]["GRCh38 Confidently mapped read pairs"] < 80:
            Warnings = Warnings + "GRCh38 Confidently mapped read pairs < 80%, "
            TotalWarnings += 1
elif args.genome == "GRCm39":
        if df.iloc[0]["GRCm39 Confidently mapped read pairs"] < 80:
            Warnings = Warnings + "GRCm39 Confidently mapped read pairs < 80%, "
            TotalWarnings += 1
    
if args.genome == None:
        if df.iloc[0]["Unmapped read pairs"] > 5:
            Warnings = Warnings + "Unmapped read pairs > 5%, "
            TotalWarnings += 1
elif args.genome == "GRCh38":
        if df.iloc[0]["GRCh38 Unmapped read pairs"] > 5:
            Warnings = Warnings + "GRCh38 Unmapped read pairs > 5%, "
            TotalWarnings += 1
elif args.genome == "GRCm39":
        if df.iloc[0]["GRCm39 Unmapped read pairs"] > 5:
            Warnings = Warnings + "GRCm39 Unmapped read pairs > 5%, "
            TotalWarnings += 1
    

if args.genome == None:
        if df.iloc[0]["Non-nuclear read pairs"] > 20:
            Warnings = Warnings + "Non-nuclear read pairs > 20%, "
            TotalWarnings += 1
elif args.genome == "GRCh38":
        if df.iloc[0]["GRCh38 Non-nuclear read pairs"] > 20:
            Warnings = Warnings + "GRCh38 Non-nuclear read pairs > 20%, "
            TotalWarnings += 1
elif args.genome == "GRCm39":
        if df.iloc[0]["GRCm39 Non-nuclear read pairs"] > 20:
            Warnings = Warnings + "GRCm39 Non-nuclear read pairs > 20%, "
            TotalWarnings += 1

    
if args.genome == None:
        if df.iloc[0]["Fragments in nucleosome-free regions"] < 40:
            MajorWarnings = MajorWarnings + "Fragments in nucleosome-free regions < 40%, "
            TotalWarnings += 1
elif args.genome == "GRCh38":
        if df.iloc[0]["GRCh38 Fragments in nucleosome-free regions"] < 40:
            MajorWarnings = MajorWarnings + "GRCh38 Fragments in nucleosome-free regions < 40%, "
            TotalWarnings += 1
elif args.genome == "GRCm39":
        if df.iloc[0]["GRCm39 Fragments in nucleosome-free regions"] < 40:
            MajorWarnings = MajorWarnings + "GRCm39 Fragments in nucleosome-free regions < 40%, "
            TotalWarnings += 1
    
    
# Targeting Metrics
if args.genome == None:
        if df.iloc[0]["Number of peaks"] < 45000:
            MajorWarnings = MajorWarnings + "Number of peaks < 45000, "
            TotalWarnings += 1
elif args.genome == "GRCh38":
        if df.iloc[0]["GRCh38 Number of peaks"] < 45000:
            MajorWarnings = MajorWarnings + "GRCh38 Number of peaks < 45000, "
            TotalWarnings += 1
elif args.genome == "GRCm39":
        if df.iloc[0]["GRCm39 Number of peaks"] < 45000:
            MajorWarnings = MajorWarnings + "GRCm39 Number of peaks < 45000, "
            TotalWarnings += 1
    
if args.genome == None:
      if 2 < df.iloc[0]["Fraction of genome in peaks"] < 20:
            MajorWarnings = MajorWarnings + "Fraction of genome in peaks < 2% or > 20%, "
            TotalWarnings += 1
elif args.genome == "GRCh38":
      if 2 < df.iloc[0]["GRCh38 Fraction of genome in peaks"] < 20:
            MajorWarnings = MajorWarnings + "GRCh38 Fraction of genome in peaks < 2% or > 20%, "
            TotalWarnings += 1
elif args.genome == "GRCm39":
      if 2 < df.iloc[0]["GRCm39 Fraction of genome in peaks"] < 20:
            MajorWarnings = MajorWarnings + "GRCm39 Fraction of genome in peaks < 2% or > 20%, "
            TotalWarnings += 1
if args.genome == None:
        if df.iloc[0]["TSS enrichment score"] < 5:
            MajorWarnings = MajorWarnings + "TSS enrichment score < 5%, "
            TotalWarnings += 1
elif args.genome == "GRCh38":
        if df.iloc[0]["GRCh38 TSS enrichment score"] < 5:
            MajorWarnings = MajorWarnings + "GRCh38 TSS enrichment score < 5%, "
            TotalWarnings += 1
elif args.genome == "GRCm39":
        if df.iloc[0]["GRCm39 TSS enrichment score"] < 5:
            MajorWarnings = MajorWarnings + "GRCm39 TSS enrichment score < 5%, "
            TotalWarnings += 1
    
if args.genome == None:
        if df.iloc[0]["Fraction of high-quality fragments overlapping TSS"] < 15:
            MajorWarnings = MajorWarnings + "Fraction of high-quality fragments overlapping TSS < 15%, "
            TotalWarnings += 1
elif args.genome == "GRCh38":
        if df.iloc[0]["GRCh38 Fraction of high-quality fragments overlapping TSS"] < 15:
            MajorWarnings = MajorWarnings + "GRCh38 Fraction of high-quality fragments overlapping TSS < 15%, "
            TotalWarnings += 1
elif args.genome == "GRCm39":
        if df.iloc[0]["GRCm39 Fraction of high-quality fragments overlapping TSS"] < 15:
            MajorWarnings = MajorWarnings + "GRCm39 Fraction of high-quality fragments overlapping TSS < 15%, "
            TotalWarnings += 1
    
if args.genome == None:
        if df.iloc[0]["Fraction of high-quality fragments overlapping peaks"] < 15:
            MajorWarnings = MajorWarnings + "Fraction of high-quality fragments overlapping peaks < 15%, "
            TotalWarnings += 1
elif args.genome == "GRCh38":
        if df.iloc[0]["GRCh38 Fraction of high-quality fragments overlapping peaks"] < 15:
            MajorWarnings = MajorWarnings + "GRCh38 Fraction of high-quality fragments overlapping peaks < 15%, "
            TotalWarnings += 1
elif args.genome == "GRCm39":
        if df.iloc[0]["GRCm39 Fraction of high-quality fragments overlapping peaks"] < 15:
            MajorWarnings = MajorWarnings + "GRCm39 Fraction of high-quality fragments overlapping peaks < 15%, "
            TotalWarnings += 1
    
  
df["Warnings"] = Warnings
df["MajorWarnings"] = MajorWarnings
df["Total Warnings"] = TotalWarnings
        
MasterDF = pandas.concat([MasterDF, df])

MasterDF.to_csv( args.outdir + "QC_Summary_CellRanger_Report.tsv", sep = "\t", index = False)
