import pandas as pd

if __name__ == '__main__':

    file = r'11_score_trio.csv'
    with open(file, 'r') as f:
        data = pd.read_csv(f)

    result = pd.DataFrame()

    groups = data.groupby(by=['Accession Number'])
    for nc, group in groups:

        max_trio_score = max(group['trio score'])

        group_trio = group[group['trio score'] == max_trio_score]

        min_mismatch = min(group_trio['min mismatch'])

        group_box = group_trio[group_trio['min mismatch'] == min_mismatch]

        #
        min_spacer = min(group_box['coserved spacer len - 13'])

        group_spacer = group_box[group_box['coserved spacer len - 13'] == min_spacer]

        result = pd.concat([result, group_spacer])

    result = result.drop_duplicates(['Accession Number', 'DoriC AC', 'Upstream Box',
                                     'seq_between_box', 'Downstream Box', 'GC-rich region', 'DnaA-trio'])
    result.to_csv('12_screen.csv', index=False)
