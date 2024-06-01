import pandas as pd


def read_smap(smap_file_path):
    headers = []
    with open(smap_file_path) as fp_read:
        for i in range(6):
            fp_read.readline()
        headers = fp_read.readline().replace('\n', '').split('\t')[1:]
    df = pd.read_csv(smap_file_path, sep='\t', skiprows=8, header=None, names=headers)
    return df


def filter_confidence(input_df, threshold=0.0):
    return input_df[input_df['Confidence'] >= threshold]


def filter_chroms(input_df, chrom1, chrom2):
    return input_df[(input_df['RefcontigID1'] == chrom1) & (input_df['RefcontigID2'] == chrom2)]


def filter_event_type(input_df, event_types):
    return input_df[input_df['Type'].isin(event_types)]



c_df = read_smap('/media/zhaoyang-new/workspace/KarSim/packaged_data_for_comparison/Bionano_solve/23Y_NF1_microdeletion_r2.1/exp_refineFinal1_merged_filter_inversions.smap')
c_df = filter_confidence(c_df)
c_df = filter_chroms(c_df, 14, 9)
# c_df = filter_event_type(c_df, ['deletion'])
print(c_df)
