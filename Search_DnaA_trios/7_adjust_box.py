import pandas as pd

file = r'6_search_upstream_box.csv'
with open(file, 'r') as f:
    data = pd.read_csv(f)

result = pd.DataFrame()

for i in data.index:
    result.at[i, 'Accession Number'] = data.at[i, 'NC']
    result.at[i, 'DoriC AC'] = data.at[i, 'ORI']

    if str(data.at[i, 'box_0']) == 'nan' and str(data.at[i, 'box2 in trio']) == 'nan':
        result.at[i, 'Upstream Box'] = data.at[i, 'box1 in trio']
        result.at[i, 'GC-rich region'] = str(data.at[i, 'GC_rich_region'])
        result.at[i, 'DnaA-trio'] = data.at[i, 'trio']

        result.at[i, 'Upstream mismatch'] = data.at[i, 'box1 mismatch']

        result.at[i, 'box motif'] = data.at[i, 'box motif']
        result.at[i, 'chromosome size'] = data.at[i, 'chromosome size']
        result.at[i, 'strand'] = data.at[i, 'strand']
        result.at[i, 'the number of trio'] = data.at[i, 'the number of trio']

    elif str(data.at[i, 'box_0']) == 'nan' and str(data.at[i, 'box2 in trio']) != 'nan':
        result.at[i, 'Upstream Box'] = data.at[i, 'box1 in trio']
        result.at[i, 'seq_between_box'] = str(
            data.at[i, 'seq between box1 and box2'])
        result.at[i, 'Downstream Box'] = data.at[i, 'box2 in trio']
        result.at[i, 'GC-rich region'] = data.at[i, 'GC_rich_region']
        result.at[i, 'DnaA-trio'] = data.at[i, 'trio']

        result.at[i, 'Upstream mismatch'] = data.at[i, 'box1 mismatch']
        result.at[i, 'Downstream mismatch'] = data.at[i, 'box2 mismatch']

        result.at[i, 'box motif'] = data.at[i, 'box motif']
        result.at[i, 'chromosome size'] = data.at[i, 'chromosome size']
        result.at[i, 'strand'] = data.at[i, 'strand']
        result.at[i, 'the number of trio'] = data.at[i, 'the number of trio']

    elif str(data.at[i, 'box_0']) != 'nan' and str(data.at[i, 'box2 in trio']) == 'nan':
        result.at[i, 'Upstream Box'] = data.at[i, 'box_0']
        result.at[i, 'seq_between_box'] = data.at[i, 'box0_seq_box1']
        result.at[i, 'Downstream Box'] = data.at[i, 'box1 in trio']
        result.at[i, 'GC-rich region'] = data.at[i, 'GC_rich_region']
        result.at[i, 'DnaA-trio'] = data.at[i, 'trio']

        result.at[i, 'Upstream mismatch'] = data.at[i, 'box0 mismatch']
        result.at[i, 'Downstream mismatch'] = data.at[i, 'box1 mismatch']

        result.at[i, 'box motif'] = data.at[i, 'box motif']
        result.at[i, 'chromosome size'] = data.at[i, 'chromosome size']
        result.at[i, 'strand'] = data.at[i, 'strand']
        result.at[i, 'the number of trio'] = data.at[i, 'the number of trio']
    else:
        print(data.at[i, 'NC'])

result.to_csv('7_adjust_box.csv', index=False)
