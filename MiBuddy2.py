import argparse
import glob
import os
import subprocess
import numpy as np
import pandas as pd
from collections import OrderedDict


# Returns chains to export from file_name

def ChainFromFilename(filename):
    chain_list = ['TRA', 'TRB', 'TRG', 'TRD', 'TCR', 'IGH', 'IGK', 'IGL', 'IG']
    chains = []
    getchain = (x for x in filename.split("_") if x in chain_list)
    for x in getchain:
        chains.append(x)
    if not chains:
        return 'ALL'
    else:
        return ','.join(chains)


# report
def report(minimal_overseq, downsample):
    estim_col_merge = ["TOTAL_READS", "#SAMPLE_ID", "TOTAL_MIGS"]
    basicstat_col_merge = ['count', 'diversity', 'mean_cdr3nt_length', 'mean_insert_size', 'mean_ndn_size']
    CdrAAprofile_cols = ["mjenergy", "kf4", "volume", "strength"]
    divers_cols = ['chao1_mean', 'observedDiversity_mean', 'normalizedShannonWienerIndex_mean']

    rename_columns = OrderedDict([('TOTAL_READS', 'Total_reads'), ('TOTAL_MIGS', 'cDNA_molecules_UMI'), \
                                  ('OVERSEQ_THRESHOLD', 'Reads_per_UMI_threshold'),
                                  ('count', 'cDNA_molecules_UMI_after_filtering'), ('diversity', 'Clonotypes'), \
                                  ('mean_cdr3nt_length', 'CDR3_length'), ('mean_insert_size', 'Added_N_nucleotides'), \
                                  ('mean_ndn_size', 'NdN')])

    basicstats = pd.read_table("vdjtools/basicstats.txt")
    estimates = pd.read_table("migec/histogram/estimates.txt")
    CdrAAProfile = pd.read_table('vdjtools/cdr3aa.stat.wt.unnorm.txt')
    diversity = pd.read_table('vdjtools/downsample_' + str(downsample) + '/diversity.strict.exact.txt')

    met_cols = [col for col in basicstats if ('label' in col) or ('sample_id' in col)]

    report = basicstats.loc[:, met_cols]
    if minimal_overseq is None:
        report = report.merge(estimates[estim_col_merge + ['OVERSEQ_THRESHOLD']], left_on="sample_id",
                              right_on="#SAMPLE_ID", \
                              how="outer").drop(['#SAMPLE_ID'], axis=1)
    else:
        report = report.merge(estimates[estim_col_merge], left_on="sample_id", right_on="#SAMPLE_ID", how="outer").drop(
            ['#SAMPLE_ID'], axis=1)
        report["OVERSEQ_THRESHOLD"] = minimal_overseq

    report = report.merge(basicstats[basicstat_col_merge + ["sample_id"]], on="sample_id", how="outer")
    report = add_property(report, CdrAAProfile, CdrAAprofile_cols)
    report["Downsample_UMI"] = downsample
    report = report.merge(diversity[divers_cols + ['sample_id']], on="sample_id", how="outer")
    report = report.round(2)

    pd.options.mode.chained_assignment = None
    for row in report.itertuples():
        if row.count < row.Downsample_UMI:
            for n in divers_cols:
                report[n][row.Index] = np.NaN

    report.rename(columns=rename_columns, inplace=True)
    report.to_csv("report_simple.txt", sep='\t', index=False, na_rep="NA")

    report_concat = pd.concat({'Metadata': report[met_cols], 'Basic_Statistics': report[list(rename_columns.values())], \
                               'CDR3_AA_physical_properties': report[CdrAAprofile_cols], \
                               'Diversity_statistics': report[["Downsample_UMI"] + divers_cols]}, axis=1).reindex(
        columns=["Metadata", 'Basic_Statistics', "CDR3_AA_physical_properties", "Diversity_statistics"], level=0)
    report_concat.index += 1

    report_concat.to_excel("report.xls", header=True, na_rep="NA", merge_cells=True)


# Returns downsample UMI value

def downsample_threshold(basicstats_df):
    a = basicstats_df["count"].quantile(q=0.2) / 2
    if basicstats_df["count"].min() > a:
        x = int(np.floor(basicstats_df["count"].min() / 100) * 100)
    else:
        x = int(np.floor(a / 100) * 100)
    if x >= 500:
        return x
    else:
        return 500


# Generating parameters for MiGec assembly
def assemble_param(minimal_overseq):
    global output_dir
    samples_overseq = {}
    with open("migec/histogram/estimates.txt") as threshold:
        for line in threshold:
            if minimal_overseq is None:
                samples_overseq[line.split()[0]] = line.split()[4]
                output_dir = "assemble"
            else:
                samples_overseq[line.split()[0]] = minimal_overseq
                output_dir = "assemble_t" + str(minimal_overseq)
    return samples_overseq, output_dir


# Creating metadata file for VDJtools
def metadata_creator():
    label_list = []
    for file in glob.glob("mixcr/*.vdjca"):
        file_label_list = []
        file_id = os.path.splitext(os.path.basename(file))[0]
        file_label_list.append(file_id + ".txt")
        file_label_list.append(file_id)
        file_label_list.extend(file_id.split("_"))
        label_list.append(file_label_list)
    maxLen = max(len(l) for l in label_list)
    metadata = pd.DataFrame(label_list)
    col_names = ["#file name", "sample_id", "label_1"]
    for i in range(3, maxLen):
        col_names.append("label_" + str(i - 1))
    metadata.columns = col_names
    metadata.to_csv("mixcr/metadata.txt", sep='\t', index=False, na_rep="NA")


# Deletes empty samples from metadata
def metadata_drop_zero_count(downsample):
    samples_id_drop = []
    basicstats = pd.read_table("vdjtools/basicstats.txt")
    metadata_downsample = pd.read_table('vdjtools/downsample_' + str(downsample) + '/metadata.txt')
    for i in basicstats['sample_id']:
        if basicstats.loc[basicstats['sample_id'] == i, 'count'].all() == 0:
            samples_id_drop.append(i)
    metadata_filter = metadata_downsample[~metadata_downsample['sample_id'].isin(samples_id_drop)]
    metadata_filter.to_csv('vdjtools/downsample_' + str(downsample) + '/metadata_filter.txt', sep='\t', index=False)


# Adds physical properties to the report dataframe
def add_property(df, CdrAAProfile_df, property_list):
    for i in property_list:
        if i not in df:
            property_table = CdrAAProfile_df[CdrAAProfile_df['property'].str.contains(i)]
            df = pd.merge(df, property_table[['sample_id', 'mean']], on='sample_id', how='left')
            df = df.rename(columns={"mean": i})
    return df


def migec_checkout(barcodesFile):
    FNULL = open(os.devnull, 'w')
    demultiplexing = subprocess.Popen(
        ['migec', 'CheckoutBatch', '-cute', '--skip-undef', barcodesFile, 'migec/checkout/'],
        stdout=FNULL, stderr=FNULL)
    demultiplexing.wait()


def migec_histogram():
    FNULL = open(os.devnull, 'w')
    hist = subprocess.Popen(['migec', 'Histogram', 'migec/checkout/', 'migec/histogram/'], stdout=FNULL, stderr=FNULL)
    hist.wait()


def migec_assemble(file_R1, file_R2, overseq, output_dir):
    FNULL = open(os.devnull, 'w')
    assemble = subprocess.Popen(['migec', '-Xmx30G', 'Assemble', '-m', overseq, '--filter-collisions', file_R1, file_R2,
                                 "migec/" + output_dir + "/"], stdout=FNULL, stderr=FNULL)
    assemble.wait()


def mixcr_align(species, file_R1, file_R2):
    print("Starting MiXCR alignment for " + os.path.splitext(os.path.basename(file_R1))[0].split("_R1")[0])
    FNULL = open(os.devnull, 'w')
    mixcr_alignment = subprocess.Popen(['mixcr', 'align', '-r', 'mixcr/alignmentReport.txt', '-f', '-s', species,
                                        file_R1, file_R2,
                                        'mixcr/' + os.path.splitext(os.path.basename(file_R1))[0].split("_R1")[
                                            0] + '.vdjca'],
                                       stdout=FNULL, stderr=FNULL)
    mixcr_alignment.wait()


def mixcr_assemble(vdjca_file, ig):
    FNULL = open(os.devnull, 'w')
    print("Starting MiXCR assemble for " + vdjca_file)
    args_mixcr_assemple = ['mixcr', 'assemble', '-r', 'mixcr/assembleReport.txt', '-f', 'mixcr/' +
                           vdjca_file + '.vdjca', 'mixcr/' + vdjca_file + '.clns']
    if ig is True or 'IG' in vdjca_file:
        args_mixcr_assemple.insert(5, '-OseparateByC=true')

    mixcr_assemble = subprocess.Popen(args_mixcr_assemple, stdout=FNULL, stderr=FNULL)
    mixcr_assemble.wait()


def mixcr_export(clns_file):
    FNULL = open(os.devnull, 'w')
    chain = ChainFromFilename(clns_file)
    print('Exporting ' + chain + ' clones for ' + clns_file)
    mixcr_export = subprocess.Popen(
        ['mixcr', 'exportClones', '-o', '--filter-stops', '-f', '-c', chain, 'mixcr/' + clns_file + '.clns',
         'mixcr/' + clns_file + '.txt'],
        stdout=FNULL, stderr=FNULL)
    mixcr_export.wait()


# Converting mixcr output for VDJTools, calc basic stats
def vdjtools_convert():
    print("Converting files to vdjtools format")
    FNULL = open(os.devnull, 'w')
    vdjtools_convert = subprocess.Popen(['vdjtools', 'Convert', '-S', 'MiXCR', '-m', 'mixcr/metadata.txt',
                                         'vdjtools/'], stdout=FNULL, stderr=FNULL)
    vdjtools_convert.wait()


def vdjtools_CalcBasicStats():
    print("Calculating basic statistics")
    FNULL = open(os.devnull, 'w')
    vdjtools_basicstats = subprocess.Popen(['vdjtools', 'CalcBasicStats', '-m', 'vdjtools/metadata.txt', 'vdjtools/'],
                                           stdout=FNULL, stderr=FNULL)
    vdjtools_basicstats.wait()


def vdjtools_CalcCdrAAProfile():
    print("Calculating CDR AA physical properties")
    FNULL = open(os.devnull, 'w')
    vdjtools_cdr_prop = subprocess.Popen(['vdjtools', 'CalcCdrAaStats', '-a',
                                          'strength,kf10,turn,cdr3contact,rim,alpha,beta,polarity,charge,surface,hydropathy,count,mjenergy,volume,core,disorder,kf2,kf1,kf4,kf3,kf6,kf5,kf8,kf7,kf9',
                                          '-w', '-r', 'cdr3-center-5', '-m', 'vdjtools/metadata.txt', 'vdjtools/'],
                                         stdout=FNULL, stderr=FNULL)
    vdjtools_cdr_prop.wait()


def vdjtools_DownSample(downsample):
    print("Downsampling data to " + str(downsample) + ' events')
    FNULL = open(os.devnull, 'w')
    vdjtools_downsample = subprocess.Popen(
        ['vdjtools', 'DownSample', '-x', str(downsample), '-m', 'vdjtools/metadata.txt',
         'vdjtools/downsample_' + str(downsample) + '/'], stdout=FNULL, stderr=FNULL)
    vdjtools_downsample.wait()


def vdjtools_CalcDiversityStats(downsample):
    print("Calculating diversity statistics for downsampled data")
    FNULL = open(os.devnull, 'w')
    vdjtools_diversity = subprocess.Popen(
        ['vdjtools', 'CalcDiversityStats', '-m', 'vdjtools/downsample_' + str(downsample) + '/metadata_filter.txt',
         'vdjtools/downsample_' + str(downsample) + '/'], stdout=FNULL, stderr=FNULL)
    vdjtools_diversity.wait()


def pipeline(barcodesFile, species, minimal_overseq, ig):
    print("\033[1;36;40mMiBuddy will take care of your data\033[0m")
    print("Starting demultiplexing")
    migec_checkout(barcodesFile)
    print("Demultiplexing is complete")
    print("Collecting MIG statistics")
    migec_histogram()
    print("MIG statistics has been calculated")
    samples_overseq = assemble_param(minimal_overseq)[0]
    assemble_path = assemble_param(minimal_overseq)[1]
    for file in glob.glob("migec/checkout/*_R1.fastq.gz"):
        filename = os.path.splitext(os.path.basename(file))[0].split("_R1")[0]
        if filename in samples_overseq.keys():
            print("Assembling MIGs for {0}. Minimal number of reads per MIG: {1}".format(filename, str(
                samples_overseq[filename])))
            file_1_path = "migec/checkout/" + filename + "_R1" + ".fastq.gz"
            file_2_path = "migec/checkout/" + filename + "_R2" + ".fastq.gz"
            overseq = samples_overseq[filename]
            migec_assemble(file_1_path, file_2_path, str(overseq), assemble_path)
            mixcr_align(species, glob.glob("migec/{0}/{1}_R1*.fastq".format(assemble_path, filename))[0],
                        glob.glob("migec/{0}/{1}_R2*.fastq".format(assemble_path, filename))[0])
            mixcr_assemble(filename, ig)
            mixcr_export(filename)

    print("Creating metadata file")
    metadata_creator()
    vdjtools_convert()
    vdjtools_CalcBasicStats()
    vdjtools_CalcCdrAAProfile()
    downsample = downsample_threshold(pd.read_table("vdjtools/basicstats.txt"))
    vdjtools_DownSample(downsample)
    metadata_drop_zero_count(downsample)
    vdjtools_CalcDiversityStats(downsample)
    print("Generating a report file")
    report(minimal_overseq, downsample)


def main(args):
    dirs = ["migec", "mixcr", "vdjtools"]
    for item in dirs:
        if not os.path.exists(item):
            os.makedirs(item)
    pipeline(args.file_with_barcodes, args.s, args.overseq, args.ig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("file_with_barcodes", help="Specify barcodes file")
    parser.add_argument("-s", help="Specify species: mmu for Mus musculus, hsa - Homo sapiens")
    parser.add_argument("-ig", help="Separate IG clones by isotypes",
                        action='store_true')
    parser.add_argument("--overseq", "-minimal_overseq", type=int, default=None,
                        help="Force minimal overseq value for all samples")
    args = parser.parse_args()
    main(args)
