import pandas as pd

file = r'7_adjust_box.csv'
with open(file, 'r') as f:
    data = pd.read_csv(f)

for i in data.index:
    # print(i)
    if str(data.at[i, 'seq_between_box']) == 'nan':
        seq_between_box = ''
    else:
        seq_between_box = data.at[i, 'seq_between_box']

    if str(data.at[i, 'Downstream Box']) == 'nan':
        box2 = ''
    else:
        box2 = data.at[i, 'Downstream Box']

    if str(data.at[i, 'GC-rich region']) == 'nan':
        gc_rich = ''
    else:
        gc_rich = data.at[i, 'GC-rich region']

    spacer = seq_between_box + box2 + gc_rich

    data.at[i, 'up_seq_trio'] = spacer
    data.at[i, 'up_seq_trio_len'] = len(spacer)

    if len(data.at[i, 'Upstream Box']) == 8:
        if str(data.at[i, 'Downstream Box']) == 'nan':
            print(data.at[i, 'NC'])
        else:
            data.at[i, 'Upstream Box'] = data.at[i, 'Upstream Box'] + \
                data.at[i, 'Downstream Box'][0]
            data.at[i, 'Downstream Box'] = data.at[i, 'Downstream Box'][1:]
            data.at[i, 'up_seq_trio'] = spacer[1:]
            data.at[i, 'up_seq_trio_len'] = len(spacer)-1

    data.at[i, 'down_seq_trio'] = gc_rich
    data.at[i, 'down_seq_trio_len'] = len(gc_rich)

data.to_csv('8_calculate_spacer_len.csv', index=False)
