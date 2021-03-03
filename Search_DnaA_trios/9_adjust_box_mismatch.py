import pandas as pd

file = r'8_calculate_spacer_len.csv'
with open(file, 'r') as f:
    data = pd.read_csv(f)

for i in data.index:
    box_motif = data.at[i, 'box motif']

    if '[' in box_motif:
        data.at[i, 'Upstream mismatch'] = data.at[i, 'Upstream mismatch'] + 0.5
        data.at[i, 'Downstream mismatch'] = data.at[i,
                                                    'Downstream mismatch'] + 0.5

data.to_csv('9_adjust_box_mismatch.csv', index=False)
