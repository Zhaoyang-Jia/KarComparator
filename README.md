# KarComparator

## Pipeline: Running OMKar and comparing
All code will be run under path `KarComparator/omkar_analyses_pipeline/`

### Modify Environmental Variables:
*Please see beginning of `run_and_analyze_omkar.sh`*

**data_dir**: absolute path to DIR containing Bionano's outputs, should be individual folders in this DIR, 
where each folder contains the cnv, smap, rcmap, and xmap files of a SINGLE case.

**exe_dir**: absolute path to OMKar's DIR

### Running the Pipeline

`bash run_and_analyze_omkar.sh [version]`

For example, `bash run_and_analyze_omkar.sh b14`.

### Outputs
Each run's output will be in `KarComparator/omkar_analyses_pipeline/builds/<version>`
1. A commit version recording in `omkar_versions.txt`, associating the build name with the current OMKar commit code 
2. Intermediate Files
   1. OMKar's total output, including debug logs
   2. Dependent Cluster files
3. Summary statistics: `analyses_summary/summary_statistics.txt`
4. CSV file of all cases with different number of paths: `analyses_summary/cases_with_path_diff.csv`
5. CSV file of all cases with missed SVs: `analyses_summary/cases_with_missed_SV.csv`
6. CSV file of all statistics: `analyses_summary/whole_analyses_table.csv`

| column name                                          | explanation                                                   |
|------------------------------------------------------|---------------------------------------------------------------|
| file_name                                            | simulation case name                                          |
| cluster                                              | the number of the cluster corresponding to the cluster file   |
| n_origin_chr                                         | number of origin chromosomes for the segments in this cluster |
| origin_chr                                           | the origin chromosomes for the segments in this cluster       |
| n_path_karsim/n_path_omkar                           | number of paths in karsim/omkar                               |
| event_chr                                            | origin chromosomes that was recorded in karsim to have events |
| histories                                            | all events recorded in karsim that happened on this cluster   |
| n_path_diff                                          | #OMKar_Path - #KarSim_Path                                    |
| karsim_missed_transition(bad naming)                 | SV edges that are missed by OMKar                             |
| n_karsim_missed_transition/n_omkar_missed_transition | number of SV edges missed by OMKar/added by OMKar             |
| n_preILP_similar_SV                                  | number of missed SV edges that are present in preILP          |
| karsim_missed_transition_event_type                  | missed SV edges' corresponding event recorded in Karsim       |


## Pipeline: Compare cases between two runs

### Running the Pipeline

`python compare_two_runs.py [version1] [version2]`

For example, `python compare_two_runs.py b12 b14`

### Outputs
Each run's output will be in `KarComparator/omkar_analyses_pipeline/comparisons`

1. CSV of cases where "n_different_path" changed between the versions: `<version1>vs<version2>_path_diffs.csv`
   1. "n_diff_path" = #OMKar_Path - #KarSim_Path
   2. For example, if n_path_diff_df1 = 1 and n_path_diff_df2 = 2, 
   3. this means we used to report 1 additional path in version1, 
   4. but we are now reporting 2 additional paths in version2
2. CSV of cases where number of SV missed changed between the versions: `<version1>vs<version2>_path_diffs.csv`
   1. "SV_missed_delta" = #Version2_SV_missed - #Version1_SV_missed

## (Pending) Visualization: Make preILP/postILP graph fior OMKar
