import Analyses_UTILS as AU
import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--data_dir', type=str, help='path to DIR of dependent cluster files')
parser.add_argument('--output_dir', type=str, help='path to DIR of analyses outputs')
parser.add_argument('--omkar_log_dir', type=str, help='path to DIR of OMKar total output')
args = parser.parse_args()

data_dir = args.data_dir
output_dir = args.output_dir
omkar_log_dir = args.omkar_log_dir
AU.data_folder = data_dir
AU.karsim_file_prefix = '../new_data_files/KarSimulator/'
AU.forbidden_region_file = '../Metadata/acrocentric_telo_cen.bed'
AU.omkar_log_dir = omkar_log_dir
AU.karsim_history_edges_folder = '../packaged_data/Karsimulator_history_intermediate_terminal_labeled/'
cn_bin_file = '../Metadata/cn_bins_50kbp.txt'

df = AU.prep_df()
df = AU.process_comparison(df)
df[['preILP_similar_SV', 'n_preILP_similar_SV']] \
    = df.apply(lambda row: pd.Series(AU.iterative_check_missed_SV_in_preILP(row['file_name'], row['karsim_missed_transition'])), axis=1)
df = AU.label_missed_SV_edges(df)

# CN
print('CN started')
df[['karsim_CN', 'omkar_CN']] \
    = df.apply(lambda row: pd.Series(AU.iterative_get_cn_with_bins(row, cn_file_name=cn_bin_file)), axis=1)
print('CN conglomerate started')
karsim_sum_CN = df.groupby('file_name')['karsim_CN'].apply(AU.vector_addition)
karsim_sum_CN_dict = karsim_sum_CN.to_dict()
omkar_sum_CN = df.groupby('file_name')['omkar_CN'].apply(AU.vector_addition)
omkar_sum_CN_dict = omkar_sum_CN.to_dict()
# print('cos sim started')
# case_cos_sim, avg_case_cos_sim = AU.cos_sim_by_case(karsim_sum_CN_dict, omkar_sum_CN_dict)
# stacked_cos_sim = AU.cos_sim_stacked(karsim_sum_CN_dict, omkar_sum_CN_dict)
print('CN Jaccard similairty started')
case_jaccard_sim = AU.cn_jaccard_sim_by_case(karsim_sum_CN_dict, omkar_sum_CN_dict)
altered_bin_union_cos_sim = AU.cos_sim_with_jaccard_union_by_case(karsim_sum_CN_dict, omkar_sum_CN_dict)
AU.plot_cn_jaccared_sim_by_case(case_jaccard_sim, output_dir + 'cn_jaccard_plot.png')
AU.output_cn_dict_by_case(case_jaccard_sim, output_dir + 'cn_jaccard_scores.txt')
AU.output_cn_dict_by_case(altered_bin_union_cos_sim, output_dir + 'cn_cos_scores.txt')
# sanity check
bp_CNV = df.groupby('file_name')['CNV_missed'].sum()
bp_CNV_dict = bp_CNV.to_dict()
for c_file_name, bp_cnv_value in bp_CNV_dict.items():
    c_case_cos_sim = altered_bin_union_cos_sim[c_file_name]
    c_case_jaccard_sim = case_jaccard_sim[c_file_name]
    karsim_cn_vector = np.array(karsim_sum_CN_dict[c_file_name])
    omkar_cn_vector = np.array(omkar_sum_CN_dict[c_file_name])
    sum_bin_abs_diff = np.sum(np.abs(omkar_cn_vector - karsim_cn_vector))
    # print('{}, cos_sim: {}, bp_CNV: {}, bin_diff: {}'.format(c_file_name, case_cos_sim[c_file_name], bp_cnv_value, sum_bin_abs_diff))
    print('{}, jaccard_sim: {}, cos_sim: {}, bp_CNV: {}, bin_diff: {}'.format(c_file_name, c_case_jaccard_sim, c_case_cos_sim, bp_cnv_value, sum_bin_abs_diff))
# print('avg case cos sim: {}'.format(avg_case_cos_sim))
# print('stacked cos sim: {}'.format(stacked_cos_sim))

n_path_diff = (df['n_path_diff'] != 0).sum()
n_initial_SV = df['n_initial_SVs'].sum()
SV_missed = df['n_karsim_missed_transition'].sum()
SV_added = df['n_omkar_missed_transition'].sum()
n_preILP_similar_SV = df['n_preILP_similar_SV'].sum()
total_dummies = df['n_dummies'].sum()
total_leftover_dummies = df['n_leftover_dummies'].sum()
n_seg_changed = df['n_seg_changed'].sum()
n_leftover_segchange = df['n_leftover_segchange'].sum()
total_significant_dummies = df['n_significant_dummies'].sum()
total_leftover_significant_dummies = df['n_leftover_significant_dummies'].sum()
signif_seg_changed = df['n_significant_seg_change'].sum()
leftover_signif_seg_changed = df['n_leftover_significant_seg_change'].sum()

correct_sig_dummies = total_significant_dummies - total_leftover_significant_dummies
correct_cn_change = signif_seg_changed - leftover_signif_seg_changed
if total_significant_dummies == 0:
    percent_sig_dummies_leftover = 0
else:
    percent_sig_dummies_leftover = correct_sig_dummies / total_significant_dummies * 100
if signif_seg_changed == 0:
    percent_sig_cn_leftover = 0
else:
    percent_sig_cn_leftover = correct_cn_change / signif_seg_changed * 100


jaccard_score = (n_initial_SV - SV_missed) / (n_initial_SV + SV_added)
precision = (n_initial_SV - SV_missed) / (n_initial_SV - SV_missed + SV_added)
recall = (n_initial_SV - SV_missed) / n_initial_SV

with open(output_dir + 'summary_statistics.txt', 'w') as fp_write:
    fp_write.write('number of cluster with number of path being different: {}\n'.format(n_path_diff))
    fp_write.write('number of initial SV: {}\n'.format(n_initial_SV))
    fp_write.write('number of SV missed: {}\n'.format(SV_missed))
    fp_write.write('number of SV added: {}\n'.format(SV_added))
    fp_write.write('number of SV missed that was found in preILP graph: {}\n'.format(n_preILP_similar_SV))
    fp_write.write('Jaccard Score: {}\n'.format(jaccard_score))
    fp_write.write('precision: {}\n'.format(precision))
    fp_write.write('recall: {}\n'.format(recall))
    fp_write.write('number of significant dummies (total): {} ({})\n'.format(total_significant_dummies, total_dummies))
    fp_write.write('number of significant dummies correct: {} ({}%)\n'.format(correct_sig_dummies, percent_sig_dummies_leftover))
    fp_write.write('number of significant CN change (total): {} ({})\n'.format(signif_seg_changed, n_seg_changed))
    fp_write.write('number of significant CN change correct: {} ({}%)\n'.format(correct_cn_change, percent_sig_cn_leftover))

path_diff_df = df[df['n_path_diff'] != 0]
path_diff_df = path_diff_df.sort_values(by=['n_path_diff'])
path_diff_df.to_csv(output_dir + 'cases_with_path_diff.csv', index=False)

sv_missed_df = df[df['n_karsim_missed_transition'] != 0]
sv_missed_df = sv_missed_df.sort_values(by=['n_karsim_missed_transition'])
sv_missed_df.to_csv(output_dir + 'cases_with_missed_SV.csv', index=False)

df.to_csv(output_dir + 'whole_analyses_table.csv', index=False)


## analyze by case
sums_by_file = df.groupby('file_name').agg({'n_initial_SVs': 'sum', 'SV_missed': 'sum', 'histories': AU.sum_history_dicts})
sums_by_file.to_csv(output_dir + 'summed_case_statistics.csv')

