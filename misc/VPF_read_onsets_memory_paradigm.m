clear;clc;

onset_table = readtable('/media/pfaffenrot/Elements/postdoc/projects/hippocampus_VASO/derivatives/pipeline/7568/ses-01/func/sub09_ses01_run3_onset.csv');
TR_vol = 2.2154;
onsets_before_press = table2array(onset_table(:,3));

%take into account the fact that we have disrecarded the first 3 volumes as dummies
offset = onsets_before_press(1)-2*TR_vol;

onsets_before_press = onsets_before_press-offset;

onsets_after_press = table2array(onset_table(:,4));
onsets_after_press = onsets_after_press-offset;

duration_before_press = table2array(onset_table(:,5));
duration_after_press = table2array(onset_table(:,6));



memory_task_idx = strcmp(onset_table.('condition'),'AM');
math_task_idx = strcmp(onset_table.('condition'),'MA');

onsets_before_press_memory = onsets_before_press(memory_task_idx);
onsets_after_press_memory = onsets_after_press(memory_task_idx);
nandex = find(isnan(onsets_after_press_memory));

onsets_before_press_math = onsets_before_press(math_task_idx);
onsets_after_press_math = onsets_after_press(math_task_idx);


duration_before_press_memory = duration_before_press(memory_task_idx);
duration_after_press_memory = duration_after_press(memory_task_idx);
duration_before_press_math = duration_before_press(math_task_idx);
duration_after_press_math = duration_after_press(math_task_idx);



onsets_after_press_memory(nandex) = [];
duration_before_press_memory(nandex)=17.6;
duration_after_press_memory(nandex) = [];
