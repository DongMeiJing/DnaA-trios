import pandas as pd

file1 = r'1_search_standard_box_spacer_0_16_greedy.csv'
file2 = r'2_search_specific_box_spacer_0_16_greedy.csv'
file3 = r'3_search_Epsilonproteobacteria_box_spacer_0_16_greedy.csv'

with open(file1, 'r') as f1:
    data1 = pd.read_csv(f1)

with open(file2, 'r') as f2:
    data2 = pd.read_csv(f2)

with open(file3, 'r') as f3:
    data3 = pd.read_csv(f3)

result = pd.concat([data1, data2, data3])


result = result.sort_values(by=['box1 mismatch'], axis=0, ascending=True)

result = result.drop_duplicates(
    subset=['NC', 'ORI', 'box_gc_trio_seq', 'the start of trio', 'the end of trio'])

result.to_csv('4_concat_delete_repeat.csv', index=False)
