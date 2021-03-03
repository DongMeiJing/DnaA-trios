import pandas as pd

if __name__ == '__main__':
    file = r'9_adjust_box_mismatch.csv'
    with open(file, 'r') as f:
        data = pd.read_csv(f)

    for i in data.index:
        up_mis = data.at[i, 'Upstream mismatch']
        down_mis = data.at[i, 'Downstream mismatch']

        if str(down_mis) == 'nan' or up_mis <= down_mis:
            data.at[i, 'min mismatch'] = up_mis
            data.at[i, 'coserved spacer'] = data.at[i, 'up_seq_trio']
            data.at[i, 'coserved spacer len'] = data.at[i, 'up_seq_trio_len']
            data.at[i, 'coserved spacer len - 13'] = abs(
                data.at[i, 'up_seq_trio_len'] - 13)
        elif down_mis < up_mis:
            data.at[i, 'min mismatch'] = down_mis
            data.at[i, 'coserved spacer'] = data.at[i, 'down_seq_trio']
            data.at[i, 'coserved spacer len'] = data.at[i, 'down_seq_trio_len']
            data.at[i, 'coserved spacer len - 13'] = abs(
                data.at[i, 'down_seq_trio_len']-13)

    data.to_csv('10_minmismatch_newspacer.csv', index=False)
