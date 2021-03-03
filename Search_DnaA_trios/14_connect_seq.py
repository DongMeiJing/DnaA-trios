import pandas as pd
import os

file = r'12_screen.csv'
with open(file, 'r') as f:
    data = pd.read_csv(f)

for i in data.index:
    box1 = data.at[i, 'Upstream Box']
    if str(data.at[i, 'up_seq_trio']) == 'nan':
        spacer = ''
    else:
        spacer = data.at[i, 'up_seq_trio']
    trio = data.at[i, 'DnaA-trio']

    all_seq = box1 + spacer + trio

    data.at[i, 'predicted sequence'] = all_seq

data.rename(
    columns={'predicted sequence': 'Predicted BUS sequences'}, inplace=True)
data.rename(columns={'Upstream Box': 'Upstream DnaA-box'}, inplace=True)
data.rename(
    columns={'seq_between_box': 'Sequence between two DnaA-boxes'}, inplace=True)
data.rename(columns={'Downstream Box': 'Downstream DnaA-box'}, inplace=True)
data.rename(columns={'trio score': 'The score of DnaA-trios'}, inplace=True)
data.rename(
    columns={'the number of trio': 'The number of DnaA-trios'}, inplace=True)
data.rename(
    columns={'box motif': 'The standard motif of DnaA-box'}, inplace=True)
data.rename(
    columns={'Upstream mismatch': 'The mismatch of upstream DnaA-box'}, inplace=True)
data.rename(columns={
            'Downstream mismatch': 'The mismatch of downstream DnaA-box'}, inplace=True)
data.rename(
    columns={'GC-rich region': 'Gap region (previous GC-rich region)'}, inplace=True)
data.rename(
    columns={'coserved spacer len': 'The length of spacer'}, inplace=True)


data = data[['Accession Number', 'DoriC AC', 'Predicted BUS sequences', 'Upstream DnaA-box', 'Downstream DnaA-box', 'DnaA-trio', 'The score of DnaA-trios', 'The number of DnaA-trios',
             'The standard motif of DnaA-box', 'The mismatch of upstream DnaA-box', 'The mismatch of downstream DnaA-box', 'Gap region (previous GC-rich region)', 'The length of spacer']]

os.remove('1_search_standard_box_spacer_0_16_greedy.csv')
os.remove('2_search_specific_box_spacer_0_16_greedy.csv')
os.remove('3_search_Epsilonproteobacteria_box_spacer_0_16_greedy.csv')
os.remove('4_concat_delete_repeat.csv')
os.remove('5_box1_box2.csv')
os.remove('6_search_upstream_box.csv')
os.remove('7_adjust_box.csv')
os.remove('8_calculate_spacer_len.csv')
os.remove('9_adjust_box_mismatch.csv')
os.remove('10_minmismatch_newspacer.csv')
os.remove('11_score_trio.csv')
os.remove('12_screen.csv')


data.to_csv('DnaA_trios.csv', index=False)
