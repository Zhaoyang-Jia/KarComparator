import Analyses_UTILS as AU
import argparse
import pandas as pd

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

df = AU.prep_df()
df = AU.process_comparison(df)
df[['preILP_similar_SV', 'n_preILP_similar_SV']] \
    = df.apply(lambda row: pd.Series(AU.iterative_check_missed_SV_in_preILP(row['file_name'], row['karsim_missed_transition'])), axis=1)
df = AU.label_missed_SV_edges(df)

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

with open(output_dir + 'summary_statistics.txt', 'w') as fp_write:
    fp_write.write('number of cluster with number of path being different: {}\n'.format(n_path_diff))
    fp_write.write('number of initial SV: {}\n'.format(n_initial_SV))
    fp_write.write('number of SV missed: {}\n'.format(SV_missed))
    fp_write.write('number of SV added: {}\n'.format(SV_added))
    fp_write.write('number of SV missed that was found in preILP graph: {}\n'.format(n_preILP_similar_SV))
    fp_write.write('Jaccard Score: {}\n'.format(jaccard_score))
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
