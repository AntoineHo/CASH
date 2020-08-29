#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Small GATK alignment and variant calling pipeline using python"""

#import errno
import os
import sys
#import shutil
import datetime
from time import localtime, strftime
import argparse
#import binascii
#import gzip
#import subprocess
#from multiprocessing import Pool, TimeoutError

import pandas as pd
import numpy as np

from Bio import SeqIO # Need BIOPYTHON SEQ/IO
from pysam import VariantFile # Need Pysam

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
#from matplotlib.lines import Line2D


#     _____   _                __    __   _                           _                  _   _
#    / ____| | |              / _|  / _| | |                         (_)                | | (_)
#   | (___   | |__    _   _  | |_  | |_  | |   ___   _ __     _ __    _   _ __     ___  | |  _   _ __     ___
#    \___ \  | '_ \  | | | | |  _| |  _| | |  / _ \ | '__|   | '_ \  | | | '_ \   / _ \ | | | | | '_ \   / _ \
#    ____) | | | | | | |_| | | |   | |   | | |  __/ | |      | |_) | | | | |_) | |  __/ | | | | | | | | |  __/
#   |_____/  |_| |_|  \__,_| |_|   |_|   |_|  \___| |_|      | .__/  |_| | .__/   \___| |_| |_| |_| |_|  \___|
#                                                            | |         | |
#                                                            |_|         |_|

def shuffler(args) :
    """Reads VCF, Reference and BED files, then shuffles sites"""
    # Get reference, fai, dict and intervals
    ref = check_files(args.FASTA)[0]
    #refindex = check_files([ref + ".fai"])[0] # UNUSED
    VCF = check_files(args.VCF)[0]
    BED = check_files(args.BED)[0]
    OUT = os.path.abspath(args.OUT[0])

    # Get other arguments
    dc_args = {"sample":args.sample[0],
               "min_cov":args.low_coverage[0],
               "max_cov":args.high_coverage[0],
               "min_qual":args.min_qual[0],
               "plot":args.plot,
               "size":args.contig_size[0]
               }

    print("# shuffle.py shuffle")
    print("Reference genome:\t{}\n".format(ref))
    print("Variants:\t{}\n".format(VCF))
    print("Regions:\t{}\n".format(BED))
    print("Other arguments: " + str(dc_args))
    print("===============================================================================\n")

    if os.path.isdir(OUT) :
        print("WARNING: Output directory already exists! {}".format(OUT))
    else :
        os.makedirs(OUT) # Create directory following path

    log("Reading fasta records...")
    ref_dict = {} # Read sequence file .fa
    for fasta_record in SeqIO.parse(open(ref, "r"), "fasta") :
        ref_dict[fasta_record.id] = len(fasta_record.seq)
    log("Done!")

    log("Reading regions bed file...")
    regions = read_bed(BED)
    print(regions)
    log("Done!")

    log("Reading VCF records...")
    variants = read_vcf(VCF, dc_args["sample"]) # Read variants file .vcf
    for chrom in variants["CHROM"].unique() :
        print("Found: {} variants on {}".format(len(variants.loc[variants["CHROM"] == chrom]), chrom))
    log("Done!")

    log("Filtering VCF records...")
    variants = filter_variants(variants, dc_args) # Read variants file .vcf
    for chrom in variants["CHROM"].unique() :
        print("Filtered: {} variants on {}".format(len(variants.loc[variants["CHROM"] == chrom]), chrom))
    print(variants)
    log("Done!")

    log("Shuffling per chromosome...")
    df = shuffle_by_chromosome(variants, regions, OUT, ref_dict, dc_args["sample"])
    log("Done!")

    if dc_args["plot"] :
        log("Outputting graphs to {}...".format(os.path.join(OUT, dc_args["sample"] + "_<CHROM>.png")))
        plot_df(df, ref_dict, OUT, sample=dc_args["sample"], minimum_contig_size=dc_args["size"])
        log("Done!")

    log("Finished!")
    sys.exit(0)


def shuffle_by_chromosome(variants, regions_df, output, ref_dict, sample) :
    """Create permutations independantly in regions and in the rest of the genome while conserving the same number of SNPs"""
    final_dataframes = {} # Storing the dataframes per chromosome

    # For each chromosome found in variants
    for chrom in variants["CHROM"].unique() :
        # Extract a 1d array containing positions of variants along a chromosome
        positions = variants.loc[variants["CHROM"] == chrom]["POS"].values

        # Create a 1d array of 0s of length = length current chromosome
        presnp = np.zeros(ref_dict[chrom])
        prereg = np.zeros(ref_dict[chrom])

        # For proofing that all the chromosome is reconstructed properly in the end with the SNPs positions shuffled independantly in REG and in the rest of the genome
        print("\t- Contig: {}".format(chrom))
        print("\t- Total length (bp): {}".format(len(presnp)))

        # Modify all positions in the SNP vector where a variant is found in the VCF so it becomes 1 instead of 0
        presnp[positions] = 1

        # Modify the prereg vector so that regions become 1 instead of 0 and not regions stay 0
        regions = regions_df.loc[regions_df["CHROM"] == chrom]
        for n, region in regions.iterrows() :
            prereg[region["START"]:region["END"]] = 1

        # Create a dataframe with POS SNP and REGIONS : False indicates nothing at this site and True indicates a SNP or if in REG
        df = pd.DataFrame.from_dict({"SNP" : presnp, "REGIONS":prereg})
        df["SNP"] = df["SNP"].astype(bool) # Convert 0.0 and 1.0 into False and True
        df["REGIONS"] = df["REGIONS"].astype(bool) # Convert 0.0 and 1.0 into False and True
        df = df.reset_index(drop=False).rename(columns={"index":"POS"}) # Reset index and rename it POS

        # Extract regions NOT in REG
        outREG = df.loc[df["REGIONS"] == False]
        # Shuffle values in SNP column of outREG while keeping index in place
        # Create a new column "shuffledSNP" containing the shuffled version
        outREG = outREG.assign(shuffledSNP=np.random.permutation(outREG["SNP"].values))
        print("\t- Length outside regions (bp): {}".format(len(outREG)))

        # Shuffle values in SNP column of each individual REG
        REGgroups = []
        insideREGbp = 0
        for name, group in df.groupby([(df["REGIONS"] != df["REGIONS"].shift()).cumsum()]) :
            if group.iloc[0,2] : # Get first row and 3rd column of group : if True == inside region
                outgroup = group.assign(shuffledSNP=np.random.permutation(group["SNP"].values))
                REGgroups.append(outgroup)
                insideREGbp += len(outgroup)

        print("\t- Length inside regions (bp): {}".format(insideREGbp))

        # Merge REG and outREG groups
        final_REG_shuffled_separately = pd.concat([outREG] + REGgroups).sort_index()
        print("\t- Concatenated length (bp): {}".format(len(final_REG_shuffled_separately)))
        #print(final_REG_shuffled_separately) # FOR DEBUG
        final_dataframes[chrom] = final_REG_shuffled_separately # Add df to a dictionnary containing DF
        print("")

    # Output to csv
    output_dfs = []
    for chrom, df in final_dataframes.items() :
        df = df.assign(CHROM=[chrom for i in range(len(df))])
        output_dfs.append(df)
    output_df = pd.concat(output_dfs)
    output_df.to_csv(os.path.join(output, "{}_DataFrame.csv".format(sample)))

    return output_df

#    _____    _           _     _     _
#   |  __ \  | |         | |   | |   (_)
#   | |__) | | |   ___   | |_  | |_   _   _ __     __ _
#   |  ___/  | |  / _ \  | __| | __| | | | '_ \   / _` |
#   | |      | | | (_) | | |_  | |_  | | | | | | | (_| |
#   |_|      |_|  \___/   \__|  \__| |_| |_| |_|  \__, |
#                                                  __/ |
#                                                 |___/

def plotter(args) :
    ref = check_files(args.FASTA)[0]
    CSV = check_files(args.CSV)[0]
    outdir = os.path.join(os.getcwd(), args.output[0]) # Outdir path

    # Get other arguments
    dc_args = {"width":args.width[0], "height":args.height[0], "bins":args.bins[0], "size":args.size[0],
               "color_observed":args.color_observed[0], "color_shuffled":args.color_shuffled[0],
               }

    print("# shuffle.py plot")
    print("Reference genome:\t{}\n".format(ref))
    print("DataFrame:\t{}\n".format(CSV))
    print("Output directory:\t{}\n".format(outdir))
    print("Other arguments: " + str(dc_args))
    print("===============================================================================\n")

    if os.path.isdir(outdir) :
        print("WARNING: Output directory already exists! {}".format(outdir))
    else :
        os.makedirs(outdir) # Create directory following path

    log("Reading fasta records...")
    ref_dict = {} # Read sequence file .fa
    for fasta_record in SeqIO.parse(open(ref, "r"), "fasta") :
        ref_dict[fasta_record.id] = len(fasta_record.seq)
    log("Done!")

    log("Loading .csv DataFrame...")
    all_df = pd.read_csv(CSV, sep=",", header=0, index_col=0, names=["POS", "SNP", "REGIONS", "shuffledSNP", "CHROM"])
    log("Done!")

    log("Plotting DataFrame...")
    plot_df(all_df, ref_dict, outdir, minimum_contig_size=dc_args["size"],
            width=dc_args["width"], height=dc_args["height"],
            nbins=dc_args["bins"], color_shuffled=dc_args["color_shuffled"],
            color_observed=dc_args["color_observed"])
    log("Done!")

    log("Finished!")
    sys.exit(0)


def plot_df(all_df, ref_dict, outdir, sample="unknown", minimum_contig_size=1000, width=14, height=7, nbins=501, color_shuffled="orange", color_observed="dodgerblue") :
    """Plot the distributions of heterozygous sites and the shuffled distribution from a generated DataFrame"""
    for chrom in all_df["CHROM"].unique() :
        # Get df (only 1 chrom)
        df = all_df.loc[all_df["CHROM"] == chrom]

        # Get current chromosome length
        clen = ref_dict[chrom]

        # Skip too small chromosomes
        if clen < minimum_contig_size :
            continue

        # Unshuffled SNP closest neighbors
        unshuf = df.loc[df["SNP"] != False]
        unshuf = unshuf.assign(CN=unshuf["POS"].rolling(window=3, min_periods=1).apply(get_closest_neighbor))

        # shuffled SNP closest neighbors
        shuf = df.loc[df["shuffledSNP"] != False]
        shuf = shuf.assign(CN=shuf["POS"].rolling(window=3, min_periods=1).apply(get_closest_neighbor))

        # Create bins/windows over the chromosome and attribute each line of the dataframes to a bin
        bins = np.linspace(0, clen, nbins)
        unshuf = unshuf.assign(BIN=pd.cut(unshuf["POS"], bins))
        unshuf_grouped = unshuf.groupby("BIN").agg({"POS":"count", "CN":"mean"})
        shuf = shuf.assign(BIN=pd.cut(shuf["POS"], bins))
        shuf_grouped = shuf.groupby("BIN").agg({"POS":"count", "CN":"mean"})

        # Initialize figure and axes
        fig, ax = plt.subplots(ncols=1,nrows=3,sharex=True,figsize=(width,height))

        # Plot #He site per bin/window
        x = np.arange(len(unshuf_grouped))
        y = unshuf_grouped["POS"]
        ax[0].plot(x, y, color="k", alpha=1.0, zorder=10)
        ax[0].fill_between(x, 0, y, color=color_observed, alpha=0.7, zorder=10, label="Observed")

        y = shuf_grouped["POS"]
        ax[0].plot(x, y, color="k", alpha=1.0, zorder=10)
        ax[0].fill_between(x, 0, y, color=color_shuffled, alpha=0.7, zorder=10, label="Shuffled")
        ax[0].set_ylabel("#He", fontsize=14)

        # Plot CN distance observed
        y = unshuf_grouped["CN"]
        median = y.median()
        ax[1].plot(x, y, color="k", alpha=1.0, zorder=10) # Raw
        ax[1].fill_between(x, 0, y, color=color_observed, alpha=0.7, zorder=10, label="Observed")
        ax[1].set_ylabel("Average\nclosest\ndistance", fontsize=14)
        ax[2].plot(x, y, color="k", alpha=1.0, zorder=10) # Zoomed
        ax[2].fill_between(x, 0, y, color=color_observed, alpha=0.7, zorder=10, label="Observed")
        ax[2].set_ylabel("Average\nclosest\ndistance\n(zoomed)", fontsize=14)

        # Plot CN distance shuffled
        y = shuf_grouped["CN"]
        median = y.median()
        ax[1].plot(x, y, color="k", alpha=1.0, zorder=10) # Raw
        ax[1].fill_between(x, 0, y, color=color_shuffled, alpha=0.7, zorder=10, label="Shuffled")
        ax[2].plot(x, y, color="k", alpha=1.0, zorder=10) # Zoomed
        ax[2].fill_between(x, 0, y, color=color_shuffled, alpha=0.7, zorder=10, label="Shuffled")
        ax[2].set_ylim(0,3*median)

        # Plot styling
        for axi in ax : # Background / ticks / grid / legend
            axi.xaxis.set_tick_params(labelsize=12, rotation=0)
            axi.yaxis.set_tick_params(labelsize=12)
            axi.set_facecolor('whitesmoke')
            axi.yaxis.grid(True, zorder=1)
            axi.xaxis.grid(True, zorder=1, which="major")
            axi.legend()
        ax[0].set_title(chrom, fontsize=16) # Title
        ax[0].set_xlim(-2, len(unshuf_grouped)+2) # X axis limits

        # Ticks styling
        minticks = []
        majticks = []
        for x in range(0, nbins) :
            if x % 10 == 0 and x % 50 != 0 :
                minticks.append(x)
            if x % 50 == 0 :
                majticks.append(x)
        ax[0].set_xticks(majticks, minor=False)
        ax[0].set_xticks(minticks, minor=True)

        ax[2].set_xlabel("Windows", fontsize=14) # X axis label

        # Plot output
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.05)
        fig.savefig(os.path.join(outdir, sample + "_{}".format(chrom) + ".png"))
        plt.close(fig)

#     _____
#    / ____|
#   | |        ___    _ __ ___    _ __     __ _   _ __    ___
#   | |       / _ \  | '_ ` _ \  | '_ \   / _` | | '__|  / _ \
#   | |____  | (_) | | | | | | | | |_) | | (_| | | |    |  __/
#    \_____|  \___/  |_| |_| |_| | .__/   \__,_| |_|     \___|
#                                | |
#                                |_|

def comparer(args) :
    """Compares Heterozygosity outside and inside regions of a BED file"""
    # Get reference, fai, dict and intervals
    ref = check_files(args.FASTA)[0]
    VCF = check_files(args.VCF)[0]
    BED = check_files(args.BED)[0]
    outdir = os.path.abspath(os.path.join(os.getcwd(), args.output[0]))

    # Get other arguments
    dc_args = {"color_in":args.color_in[0],
               "color_out":args.color_out[0],
               "width":args.width[0],
               "height":args.height[0],
               "size":args.contig_size[0],
               "region_size":args.region_size[0],
               "max_size":args.max_size[0],
               "max_het":args.max_het[0],
               "min_size":args.min_size[0],
               "min_het":args.min_het[0],
               "bin_number":args.bin_number[0],
               "bin_het_number":args.bin_het_number[0],
               "sample":args.sample[0],
               }

    print("# shuffle.py compare")
    print("Reference genome:\t{}\n".format(ref))
    print("Variants:\t{}\n".format(VCF))
    print("Regions:\t{}\n".format(BED))
    print("Other arguments: " + str(dc_args))
    print("===============================================================================\n")

    if os.path.isdir(outdir) :
        print("WARNING: Output directory already exists! {}".format(outdir))
    else :
        os.makedirs(outdir) # Create directory following path

    log("Reading fasta records...")
    ref_dict = {} # Read sequence file .fa
    for fasta_record in SeqIO.parse(open(ref, "r"), "fasta") :
        ref_dict[fasta_record.id] = len(fasta_record.seq)
    log("Done!")

    log("Reading regions bed file...")
    regions = read_bed(BED)
    print(regions)
    log("Done!")

    log("Reading VCF records...")
    variants = read_vcf(VCF, dc_args["sample"]) # Read variants file .vcf
    for chrom in variants["CHROM"].unique() :
        print("Found: {} variants on {}".format(len(variants.loc[variants["CHROM"] == chrom]), chrom))
    log("Done!")

    log("Comparing heterozygosity...")
    chrom_df, chrom_het_prop = compare_heterozygosity(variants, regions, ref_dict)
    # Output to csv
    output_dfs = []
    for chrom, df in chrom_df.items() :
        df = df.assign(CHROM=[chrom for i in range(len(df))])
        output_dfs.append(df)
    output_df = pd.concat(output_dfs)
    output_df.to_csv(os.path.join(outdir, "{}_Comparison_DataFrame.csv".format(dc_args["sample"])))
    log("Done!")

    log("Output results...")
    compare_result(chrom_df, chrom_het_prop, outdir, dc_args, ref_dict)
    log("Done!")

    log("Finished!")
    sys.exit(0)

def compare_plot(args) :
    """Compares Heterozygosity outside and inside regions of a BED file"""
    # Get reference, fai, dict and intervals
    FASTA = check_files(args.FASTA)[0]
    CSV = check_files(args.CSV)[0]
    outdir = os.path.abspath(os.path.join(os.getcwd(), args.output[0]))

    # Get other arguments
    dc_args = {"color_in":args.color_in[0],
               "color_out":args.color_out[0],
               "width":args.width[0],
               "height":args.height[0],
               "size":args.contig_size[0],
               "region_size":args.region_size[0],
               "max_size":args.max_size[0],
               "max_het":args.max_het[0],
               "min_size":args.min_size[0],
               "min_het":args.min_het[0],
               "bin_number":args.bin_number[0],
               "bin_het_number":args.bin_het_number[0],
               "sample":args.sample[0],
               }

    print("# shuffle.py compare")
    print("Reference genome:\t{}\n".format(FASTA))
    print("DataFrame:\t{}\n".format(CSV))
    print("Other arguments: " + str(dc_args))
    print("===============================================================================\n")

    if os.path.isdir(outdir) :
        print("WARNING: Output directory already exists! {}".format(outdir))
    else :
        os.makedirs(outdir) # Create directory following path

    log("Reading fasta records...")
    ref_dict = {} # Read sequence file .fa
    for fasta_record in SeqIO.parse(open(FASTA, "r"), "fasta") :
        ref_dict[fasta_record.id] = len(fasta_record.seq)
    log("Done!")

    log("Loading .csv DataFrame...")
    all_df = pd.read_csv(CSV, sep=",", header=0, index_col=0, names=["LENGTH", "REGIONS", "SNP", "CHROM"])
    log("Done!")

    log("Output results...")
    compare_plot_result(all_df, outdir, dc_args, ref_dict)
    log("Done!")

    log("Finished!")
    sys.exit(0)

def compare_result(chrom_df, chrom_het_prop, outdir, dc_args, ref_dict) :
    """Writes results to file and plot distributions"""

    # General proportion output
    f = open(os.path.join(outdir, "comparison_per_chrom.tsv"), "w")
    f.write("CHROM\t%Het_IN\t%Het_OUT\n")
    for ctg, props in chrom_het_prop.items() :
        if ctg == "TOTAL" :
            continue
        f.write("{}\t{}\t{}\n".format(ctg, props[0], props[1]))

    f.write("# Average global: {}% inside regions and {}% outside regions".format(chrom_het_prop["TOTAL"][0], chrom_het_prop["TOTAL"][1]))
    f.close()

    # Do Plots
    min_size = dc_args["size"]
    colin = dc_args["color_in"]
    colout = dc_args["color_out"]

    # BARGRAPH
    fig, ax = plt.subplots(figsize=(dc_args["width"], dc_args["height"]))
    p = 0 # Position index
    pos = [] # position for xticks
    lbl = [] # labels for xticks
    for chrom, tup in chrom_het_prop.items() :
        if chrom == "TOTAL" :
            pass # Avoid searching in the ref_dict
        elif ref_dict[chrom] < min_size : # In case too small contig
            continue
        ax.bar(height=tup, x=[p, p+1], color=[colin, colout], edgecolor="k", zorder=10)
        pos += [p+0.5]
        lbl += [chrom]
        p += 3
    # Legend format
    handles = [Patch(facecolor=colin, edgecolor="k", linewidth=1.0, label="Inside"),
              Patch(facecolor=colout, edgecolor="k", linewidth=1.0, label="Outside")]
    labels = ["Inside", "Outside"]
    ax.legend(handles=handles, labels=labels, loc="best", ncol=2, fontsize=13)
    # Xticks format
    ax.set_xticks(pos)
    ax.set_xticklabels(lbl)
    ax.xaxis.set_tick_params(labelsize=12, rotation=0)
    # Y axis format
    ax.set_ylabel("Average %Heterozygosity")
    ax.set_ylim(0,2.6)
    ax.yaxis.set_tick_params(labelsize=12)
    ax.yaxis.grid(True, zorder=1)
    # Background format
    ax.set_facecolor('whitesmoke')
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.05)
    fig.savefig(os.path.join(outdir, "bargraph_inside_outside_heterozygosity.png"))
    plt.close(fig)

    # HISTOGRAMS
    # Heterozygosity
    for chrom, df in chrom_df.items() :
        if chrom == "TOTAL" :
            pass # Avoid searching in the ref_dict
        elif ref_dict[chrom] < min_size : # In case too small contig
            continue
        cdf = df.assign(HET=df.apply(lambda x: (x["SNP"]/x["POS"])*100, axis="columns"))
        print(cdf)
        cdf = cdf.loc[cdf["POS"] >= dc_args["region_size"]] # Filter out regions smaller than "region_size"
        in_reg = cdf.loc[cdf["REGIONS"] == True]
        out_reg = cdf.loc[cdf["REGIONS"] == False]

        fig, ax = plt.subplots(figsize=(dc_args["width"], dc_args["height"]))
        ax.hist(out_reg["HET"], bins=dc_args["bin_het_number"],
                range=(dc_args["min_het"],dc_args["max_het"]), zorder=5,
                edgecolor="k", color=colout, alpha=0.7, label="Outside")
        ax.hist(in_reg["HET"], bins=dc_args["bin_het_number"],
                range=(dc_args["min_het"],dc_args["max_het"]), zorder=10,
                edgecolor="k", color=colin, alpha=0.7, label="Inside")
        # Labels
        ax.set_xlabel("%Heterozygosity", fontsize=14)
        ax.set_ylabel("#Regions", fontsize=14)
        ax.xaxis.set_tick_params(labelsize=12, rotation=0)
        ax.yaxis.set_tick_params(labelsize=12)
        ax.set_facecolor('whitesmoke')
        ax.yaxis.grid(True, zorder=1)
        ax.xaxis.grid(True, zorder=1, which="major")
        ax.set_title("All regions - {} - %Het".format(chrom), fontsize=15)
        ax.legend(loc="best", fontsize=13)
        fig.savefig(os.path.join(outdir, "{}_histogram_heterozygosity.png".format(chrom)))
        plt.close(fig)

        # Size
        fig, ax = plt.subplots(figsize=(dc_args["width"], dc_args["height"]))
        ax.hist(out_reg["POS"], bins=dc_args["bin_number"],
        range=(dc_args["min_size"],dc_args["max_size"]), zorder=5,
        edgecolor="k", color=colout, alpha=0.7, label="Outside")
        ax.hist(in_reg["POS"], bins=dc_args["bin_number"],
        range=(dc_args["min_size"],dc_args["max_size"]), zorder=10,
        edgecolor="k", color=colin, alpha=0.7, label="Inside")
        # Labels
        ax.set_xlabel("Size distribution", fontsize=14)
        ax.set_ylabel("#Regions", fontsize=14)
        ax.xaxis.set_tick_params(labelsize=12, rotation=0)
        ax.yaxis.set_tick_params(labelsize=12)
        ax.set_facecolor('whitesmoke')
        ax.yaxis.grid(True, zorder=1)
        ax.xaxis.grid(True, zorder=1, which="major")
        ax.set_title("All regions - {} - Size".format(chrom), fontsize=15)
        ax.legend(loc="best", fontsize=13)
        fig.savefig(os.path.join(outdir, "{}_histogram_heterozygosity.png".format(chrom)))
        plt.close(fig)

def compare_plot_result(chrom_df, outdir, dc_args, ref_dict) :
    """Writes results to file and plot distributions"""

    # Do Plots
    min_size = dc_args["size"]
    colin = dc_args["color_in"]
    colout = dc_args["color_out"]

    # HISTOGRAMS
    # Heterozygosity
    for chrom in chrom_df["CHROM"].unique() :
        if ref_dict[chrom] < min_size : # In case too small contig
            continue

        # Get right chromosome
        df = chrom_df.loc[chrom_df["CHROM"] == chrom]
        print(df)

        cdf = df.assign(HET=df.apply(lambda x: (x["SNP"]/x["LENGTH"])*100, axis="columns"))
        cdf = cdf.loc[cdf["LENGTH"] >= dc_args["region_size"]] # Filter out regions smaller than "region_size"
        in_reg = cdf.loc[cdf["REGIONS"] == True]
        out_reg = cdf.loc[cdf["REGIONS"] == False]

        fig, ax = plt.subplots(figsize=(dc_args["width"], dc_args["height"]))
        ax.hist(out_reg["HET"], bins=dc_args["bin_het_number"],
                range=(dc_args["min_het"],dc_args["max_het"]), zorder=5,
                edgecolor="k", color=colout, alpha=0.7, label="Outside")
        ax.hist(in_reg["HET"], bins=dc_args["bin_het_number"],
                range=(dc_args["min_het"],dc_args["max_het"]), zorder=10,
                edgecolor="k", color=colin, alpha=0.7, label="Inside")
        # Labels
        ax.set_xlabel("%Heterozygosity", fontsize=14)
        ax.set_ylabel("#Regions", fontsize=14)
        ax.xaxis.set_tick_params(labelsize=12, rotation=0)
        ax.yaxis.set_tick_params(labelsize=12)
        ax.set_facecolor('whitesmoke')
        ax.yaxis.grid(True, zorder=1)
        ax.xaxis.grid(True, zorder=1, which="major")
        ax.set_title("All regions - {} - %Het".format(chrom), fontsize=15)
        ax.legend(loc="best", fontsize=13)
        fig.savefig(os.path.join(outdir, "{}_histogram_heterozygosity.png".format(chrom)))
        plt.close(fig)

        # Size
        fig, ax = plt.subplots(figsize=(dc_args["width"], dc_args["height"]))
        ax.hist(out_reg["LENGTH"], bins=dc_args["bin_number"],
                range=(dc_args["min_size"],dc_args["max_size"]), zorder=5,
                edgecolor="k", color=colout, alpha=0.7, label="Outside")
        ax.hist(in_reg["LENGTH"], bins=dc_args["bin_number"],
                range=(dc_args["min_size"],dc_args["max_size"]), zorder=10,
                edgecolor="k", color=colin, alpha=0.7, label="Inside")
        # Labels
        ax.set_xlabel("Size distribution", fontsize=14)
        ax.set_ylabel("#Regions", fontsize=14)
        ax.xaxis.set_tick_params(labelsize=12, rotation=0)
        ax.yaxis.set_tick_params(labelsize=12)
        ax.set_facecolor('whitesmoke')
        ax.yaxis.grid(True, zorder=1)
        ax.xaxis.grid(True, zorder=1, which="major")
        ax.set_title("All regions - {} - Size".format(chrom), fontsize=15)
        ax.legend(loc="best", fontsize=13)
        fig.savefig(os.path.join(outdir, "{}_histogram_length.png".format(chrom)))
        plt.close(fig)

def compare_heterozygosity(variants, regions, ref_dict) :
    chrom_df = {} # Stores groups with counts of SNPs and size
    chrom_het_prop = {} # Stores het% inside and outside regions per chromosome
    total_in_het, total_in_bp, total_out_het, total_out_bp = 0, 0, 0, 0 # Create variables to store total % in and out

    for chrom in variants["CHROM"].unique() :
        # Extract a 1d array containing positions of variants along a chromosome
        positions = variants.loc[variants["CHROM"] == chrom]["POS"].values

        # Create a 1d array of 0s of length = chromosome
        presnp = np.zeros(ref_dict[chrom])
        prereg = np.zeros(ref_dict[chrom])

        # Modify all positions where a variant is found to be a 1
        presnp[positions] = 1

        # Extract the non-REGIONS sub-vector
        chrom_regions = regions.loc[regions["CHROM"] == chrom]
        for n, region in chrom_regions.iterrows() :
            prereg[region["START"]:region["END"]] = 1

        # Create a dataframe with POS SNP and REGIONS : False indicates nothing at this site and True indicates a SNP or if in REGIONS
        df = pd.DataFrame.from_dict({"SNP" : presnp, "REGIONS":prereg})
        df["SNP"] = df["SNP"].astype(bool) # Convert 0.0 and 1.0 into False and True
        df["REGIONS"] = df["REGIONS"].astype(bool) # Convert 0.0 and 1.0 into False and True
        df = df.reset_index(drop=False).rename(columns={"index":"POS"}) # Reset index and rename it POS

        bp_IN = len(df.loc[df["REGIONS"] == True]) # Get size of inside regions
        het_IN = len(df.loc[(df["REGIONS"] == True) & (df["SNP"] == True)]) # Get number of snps inside regions
        prop_IN = (het_IN / bp_IN) * 100 # Compute %He in region
        total_in_bp += bp_IN
        total_in_het += het_IN

        bp_OUT = len(df.loc[df["REGIONS"] == False]) # Get size of outside regions
        het_OUT = len(df.loc[(df["REGIONS"] == False) & (df["SNP"] == True)]) # Get number of snps outside regions
        prop_OUT = (het_OUT / bp_OUT) * 100 # Compute %He out region
        total_out_bp += bp_OUT
        total_out_het += het_OUT

        # Group by region
        gb = df.groupby([(df["REGIONS"] != df["REGIONS"].shift()).cumsum()])
        # For each group: store size of region (POS), if REGION is IN (True) or OUT (False) and the number of SNPs in region
        dfgb = gb.agg({"POS":"count", "REGIONS":"first", "SNP": lambda x: x.where(x == True).count()})
        dfgb = dfgb.rename({"POS":"SIZE"}) # Rename pos into size
        chrom_df[chrom] = dfgb # Store group region
        chrom_het_prop[chrom] = (prop_IN, prop_OUT) # Store %Het inside and outside regions

        print("Contig: {}".format(chrom))
        print("In regions: {}% Het".format(prop_IN))
        print("Out of regions: {}% Het".format(prop_OUT))
        print("")

    total_prop_in = (total_in_het / total_in_bp) * 100 # Compute %He out region
    total_prop_out = (total_out_het / total_out_bp) * 100 # Compute %He out region
    chrom_het_prop["TOTAL"] = (total_prop_in, total_prop_out) # Adds total proportions

    return chrom_df, chrom_het_prop


#    _____                         _
#   |  __ \                       (_)
#   | |__) |   __ _   _ __   ___   _   _ __     __ _
#   |  ___/   / _` | | '__| / __| | | | '_ \   / _` |
#   | |      | (_| | | |    \__ \ | | | | | | | (_| |
#   |_|       \__,_| |_|    |___/ |_| |_| |_|  \__, |
#                                               __/ |
#                                              |___/

def read_bed(bed) :
    df = pd.read_csv(bed, sep="\t", header=None, names=["CHROM", "START", "END"])
    df = df.assign(LENGTH=df[["START", "END"]].apply(lambda x: x["END"] - x["START"], axis="columns"))
    df = df.drop_duplicates(keep="first", ignore_index = True)
    return df

def read_vcf(vcf, sample) :
    """Reads a VCF and fetches relevant information to a Pandas DataFrame"""

    vcf_in = VariantFile(vcf)  # auto-detect input format
    vcf_in.subset_samples([sample])

    # Variant sites for probability computation (and later modelisation)
    VariantSites = {"CHROM":[], "POS":[], "TYPE":[], "DP":[], "MAF":[], "GT":[], "QUAL":[], "ALS":[]}

    for i, rec in enumerate(vcf_in) : # For each record in vcf
        if i % 200000 == 0 :
            print("Elapsed records: {}".format(i))

        gt = rec.samples[sample]["GT"] # Get sample GT
        if len(set(gt)) == 1 :
            continue # SKIP IF HOMOZYGOUS

        VariantSites["CHROM"].append(rec.chrom) # Add record CHROM
        VariantSites["POS"].append(rec.pos) # Add record POS

        als = [x for n, x in enumerate(rec.alleles) if n in gt] # Get record position
        VariantSites["ALS"].append(als) # Add record alleles
        VariantSites["GT"].append(gt) # Add sample GT

        if "<NON_REF>" in als : # In case undefined allele
            vtype = "U"
        elif any(len(x) > 1 for x in als) or "*" in als : # In case any is a deletion or is an insertion
            vtype = "I"
        else : # In case not non-ref and not an INDEL
            vtype = "S"

        VariantSites["TYPE"].append(vtype) # Add sample type

        try : # Add sample DP and compute Min AF based on allele reads frequencies. If AD or DP is unavailable fills with None
            dp = rec.samples[sample]["DP"]
            VariantSites["DP"].append(dp)
            try :
                min_ad = min(rec.samples[sample]["AD"])
                min_af = float(min_ad/dp)
                VariantSites["MAF"].append(min_af)
            except :
                VariantSites["MAF"].append(None)
        except :
            VariantSites["DP"].append(None)

        try : # Add record QUAL at this position
            VariantSites["QUAL"].append(rec.qual)
        except :
            VariantSites["QUAL"].append(None)

    return pd.DataFrame.from_dict(VariantSites)

#    ______   _   _   _                   _
#   |  ____| (_) | | | |                 (_)
#   | |__     _  | | | |_    ___   _ __   _   _ __     __ _
#   |  __|   | | | | | __|  / _ \ | '__| | | | '_ \   / _` |
#   | |      | | | | | |_  |  __/ | |    | | | | | | | (_| |
#   |_|      |_| |_|  \__|  \___| |_|    |_| |_| |_|  \__, |
#                                                      __/ |
#                                                     |___/

def filter_variants(variants, dc_args) :
    """Filter variants based on given arguments"""
    MINCOV = dc_args["min_cov"] if dc_args["min_cov"] is not None else 0
    MAXCOV = dc_args["max_cov"] if dc_args["max_cov"] is not None else 10000000
    MINQUAL = dc_args["min_qual"] if dc_args["min_qual"] is not None else 0
    return variants.loc[(variants["DP"] >= MINCOV) & (variants["DP"] <= MAXCOV) & (variants["QUAL"] >= MINQUAL)]

#     ____    _______   _    _   ______   _____     _____
#    / __ \  |__   __| | |  | | |  ____| |  __ \   / ____|
#   | |  | |    | |    | |__| | | |__    | |__) | | (___
#   | |  | |    | |    |  __  | |  __|   |  _  /   \___ \
#   | |__| |    | |    | |  | | | |____  | | \ \   ____) |
#    \____/     |_|    |_|  |_| |______| |_|  \_\ |_____/
#
#

def get_closest_neighbor(x) : # For applying a rolling window
    if len(x) == 1 :
        return int(0)
    elif len(x) == 2 :
        return(int(max(x) - min(x)))
    else :
        return int(min([x.iloc[1]-x.iloc[0], x.iloc[2]-x.iloc[1]]))

def log(string) :
    print("\n{}: {}".format(strftime("%Y-%m-%d %H:%M:%S", localtime()), string))

def check_dirs(dirs) :
    """Returns absolute paths and raise exception if dir does not exist"""
    absdirs = []
    for d in dirs :
        if not os.path.isdir(d) :
            raise Exception("ERROR: {} is not found!".format(d))
        else :
            absdirs.append(os.path.abspath(d))
    return absdirs

def check_files(files) :
    """Returns absolute file paths and raise exception if file does not exist"""
    absfiles = []
    for file in files :
        if not os.path.isfile(file) :
            raise Exception("ERROR: {} is not found!".format(file))
        else :
            absfiles.append(os.path.abspath(file))
    return absfiles

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return False

def list_str(v) :
    return v.split(',')

def str_to_bool(v) :
	if isinstance(v, bool):
		return v
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')



#               _____     _____   _    _   __  __   ______   _   _   _______    _____
#       /\     |  __ \   / ____| | |  | | |  \/  | |  ____| | \ | | |__   __|  / ____|
#      /  \    | |__) | | |  __  | |  | | | \  / | | |__    |  \| |    | |    | (___
#     / /\ \   |  _  /  | | |_ | | |  | | | |\/| | |  __|   | . ` |    | |     \___ \
#    / ____ \  | | \ \  | |__| | | |__| | | |  | | | |____  | |\  |    | |     ____) |
#   /_/    \_\ |_|  \_\  \_____|  \____/  |_|  |_| |______| |_| \_|    |_|    |_____/
#
#

def main() :
    """Argument parser"""
    parser = argparse.ArgumentParser(description='Shuffles SNPs and outputs homozygosity runs compared to contigs.')
    subparsers = parser.add_subparsers(required=True, dest="shuffle")

    # Apply shuffling
    shuffle = subparsers.add_parser('shuffle', help="Shuffles SNPs and output closest distance of observed and shuffled data.")
    shuffle.add_argument('VCF',nargs=1,type=str,help="<STRING> A path to a the short variants calls (vcf file).")
    shuffle.add_argument('FASTA',nargs=1,type=str,help="<STRING> A path to the reference genome (fasta file).")
    shuffle.add_argument('BED',nargs=1,type=str,help="<STRING> A path to the regions locations (sorted and merged bed file).")
    shuffle.add_argument('OUT',nargs=1,type=str,help="<STRING> An output directory path for the out files.")
    shuffle.add_argument('-s', '--sample',nargs=1,type=str,default=['ancestor'],help="<STRING> Sample to use for short variants positions. Default: %(default)s")
    shuffle.add_argument('-lc','--low-coverage',nargs=1,type=int,default=[None], required=False,help="<INT> Minimum coverage required to consider a good variant. Default: %(default)s")
    shuffle.add_argument('-hc','--high-coverage',nargs=1,type=int,default=[None], required=False,help="<INT> Maximum coverage threshold to consider a good variant. Default: %(default)s")
    shuffle.add_argument('-mq','--min-qual',nargs=1,type=int,default=[None], required=False,help="<INT> Minimum QUAL required to consider a good variant. Default: %(default)s")
    shuffle.add_argument('-p', '--plot',type=str_to_bool, nargs='?', const=True, default=True, help="Output distribution plots. Default: %(default)s")
    shuffle.add_argument('-cs','--contig-size',nargs=1,type=int,default=[500000], required=False, help="<INT> Minimum size of a contig to plot. Default: %(default)sbp")
    shuffle.set_defaults(func=shuffler)

    # Plot shuffled vs observed distributions based on shuffle command results
    plot = subparsers.add_parser('plot-dis', help="Plot shuffled vs observed distributions based on shuffle command results.")
    plot.add_argument('CSV',nargs=1,type=str,help="<STRING> A path to a the output of the shuffle command (csv file).")
    plot.add_argument('FASTA',nargs=1,type=str,help="<STRING> A path to the reference genome (fasta file).")
    plot.add_argument('-o', '--output',nargs=1,type=str,default=['plots'],help="<STRING> Directory name to output. Default: path/to/cur_dir/%(default)s")
    plot.add_argument('-co', '--color-observed',nargs=1,type=str,default=['dodgerblue'],help="<STRING> A matplotlib valid color for the observed distribution. Default: %(default)s")
    plot.add_argument('-cs', '--color-shuffled',nargs=1,type=str,default=['orange'],help="<STRING> A matplotlib valid color for the shuffled distribution. Default: %(default)s")
    plot.add_argument('-wd','--width',nargs=1,type=int,default=[14], required=False, help="<INT> Width of plot to output. Default: %(default)s")
    plot.add_argument('-hg','--height',nargs=1,type=int,default=[7], required=False, help="<INT> Height of plot to output. Default: %(default)s")
    plot.add_argument('-s','--size',nargs=1,type=int,default=[500000], required=False, help="<INT> Minimum size of a contig to plot. Default: %(default)sbp")
    plot.add_argument('-b','--bins',nargs=1,type=int,default=[501], required=False, help="<INT> Number of bins (windows) to make on each chromosome. Default: %(default)sbp")
    plot.set_defaults(func=plotter)

    # Compare He in and out regions of bedfile
    plot_comp = subparsers.add_parser('plot-comp', help="Compare heterozygosity of regions in and out of a bed file.")
    plot_comp.add_argument('CSV',nargs=1,type=str,help="<STRING> A path to a the short variants calls (vcf file).")
    plot_comp.add_argument('FASTA',nargs=1,type=str,help="<STRING> A path to the reference genome (fasta file).")
    plot_comp.add_argument('-s', '--sample',nargs=1,type=str,default=['unknown'],help="<STRING> Sample name (only for output file names). Default: %(default)s")
    plot_comp.add_argument('-o', '--output',nargs=1,type=str,default=['compare_plots'],help="<STRING> Directory name to output. Default: path/to/cur_dir/%(default)s")
    plot_comp.add_argument('-ci', '--color-in',nargs=1,type=str,default=['dodgerblue'],help="<STRING> A matplotlib valid color for the observed distribution. Default: %(default)s")
    plot_comp.add_argument('-co', '--color-out',nargs=1,type=str,default=['orange'],help="<STRING> A matplotlib valid color for the shuffled distribution. Default: %(default)s")
    plot_comp.add_argument('-wd','--width',nargs=1,type=int,default=[14], required=False, help="<INT> Width of plot to output. Default: %(default)s")
    plot_comp.add_argument('-hg','--height',nargs=1,type=int,default=[7], required=False, help="<INT> Height of plot to output. Default: %(default)s")
    plot_comp.add_argument('-cs','--contig-size',nargs=1,type=int,default=[500000], required=False, help="<INT> Minimum size of a contig to plot. Default: %(default)sbp")
    plot_comp.add_argument('-rs','--region-size',nargs=1,type=int,default=[200], required=False, help="<INT> Minimum size of a region to be considered in output plots. Default: %(default)sbp")
    plot_comp.add_argument('-MS','--max-size',nargs=1,type=int,default=[500000], required=False, help="<INT> Maximum size of a region to plot the histogram (regions in range [ms, MS]). Default: %(default)s")
    plot_comp.add_argument('-ms','--min-size',nargs=1,type=int,default=[200], required=False, help="<INT> Minimum size of a region to plot the histogram (regions in range [ms, MS]). Default: %(default)s")
    plot_comp.add_argument('-bn','--bin-number',nargs=1,type=int,default=[100], required=False, help="<INT> Number of bins in histogram in range [ms, MS]. Default: %(default)s")
    plot_comp.add_argument('-MH','--max-het',nargs=1,type=int,default=[10], required=False, help="<INT> Maximum heterozygosity range of a region to plot the histogram (in range [mh, MH]). Default: %(default)s")
    plot_comp.add_argument('-mh','--min-het',nargs=1,type=int,default=[0], required=False, help="<INT> Minimum heterozygosity range of a region to plot the histogram (in range [mh, MH]). Default: %(default)s")
    plot_comp.add_argument('-bhn','--bin-het-number',nargs=1,type=int,default=[100], required=False, help="<INT> Number of bins in histogram in range [mh, MH]. Default: %(default)s")
    plot_comp.set_defaults(func=compare_plot)

    # Compare He in and out regions of bedfile
    compare = subparsers.add_parser('compare', help="Compare heterozygosity of regions in and out of a bed file.")
    compare.add_argument('VCF',nargs=1,type=str,help="<STRING> A path to a the short variants calls (vcf file).")
    compare.add_argument('FASTA',nargs=1,type=str,help="<STRING> A path to the reference genome (fasta file).")
    compare.add_argument('BED',nargs=1,type=str,help="<STRING> A path to the regions locations (sorted and merged bed file).")
    compare.add_argument('-s', '--sample',nargs=1,type=str,default=['ancestor'],help="<STRING> Sample name to read in VCF. Default: %(default)s")
    compare.add_argument('-o', '--output',nargs=1,type=str,default=['plots'],help="<STRING> Directory name to output. Default: path/to/cur_dir/%(default)s")
    compare.add_argument('-ci', '--color-in',nargs=1,type=str,default=['dodgerblue'],help="<STRING> A matplotlib valid color for the observed distribution. Default: %(default)s")
    compare.add_argument('-co', '--color-out',nargs=1,type=str,default=['orange'],help="<STRING> A matplotlib valid color for the shuffled distribution. Default: %(default)s")
    compare.add_argument('-wd','--width',nargs=1,type=int,default=[14], required=False, help="<INT> Width of plot to output. Default: %(default)s")
    compare.add_argument('-hg','--height',nargs=1,type=int,default=[7], required=False, help="<INT> Height of plot to output. Default: %(default)s")
    compare.add_argument('-cs','--contig-size',nargs=1,type=int,default=[500000], required=False, help="<INT> Minimum size of a contig to plot. Default: %(default)sbp")
    compare.add_argument('-rs','--region-size',nargs=1,type=int,default=[200], required=False, help="<INT> Minimum size of a region to be considered in output plots. Default: %(default)sbp")
    compare.add_argument('-MS','--max-size',nargs=1,type=int,default=[500000], required=False, help="<INT> Maximum size of a region to plot the histogram (regions in range [ms, MS]). Default: %(default)s")
    compare.add_argument('-ms','--min-size',nargs=1,type=int,default=[200], required=False, help="<INT> Minimum size of a region to plot the histogram (regions in range [ms, MS]). Default: %(default)s")
    compare.add_argument('-bn','--bin-number',nargs=1,type=int,default=[100], required=False, help="<INT> Number of bins in histogram in range [ms, MS]. Default: %(default)s")
    compare.add_argument('-MH','--max-het',nargs=1,type=int,default=[10], required=False, help="<INT> Maximum heterozygosity range of a region to plot the histogram (in range [mh, MH]). Default: %(default)s")
    compare.add_argument('-mh','--min-het',nargs=1,type=int,default=[0], required=False, help="<INT> Minimum heterozygosity range of a region to plot the histogram (in range [mh, MH]). Default: %(default)s")
    compare.add_argument('-bhn','--bin-het-number',nargs=1,type=int,default=[100], required=False, help="<INT> Number of bins in histogram in range [mh, MH]. Default: %(default)s")
    compare.set_defaults(func=comparer)

    args = parser.parse_args()
    args.func(args)

    sys.exit(0)

if __name__ == '__main__':
	main()
