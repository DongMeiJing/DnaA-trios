import pandas as pd
import re


def get_trio_score(trio_list):
    GAT_list = ['GAT']
    AAT_list = ['AAT']
    GAX = re.compile('.AT|GA.')
    XAX = re.compile('.A.')
    trio_list_new = []
    for each_trio in trio_list:
        trio_list_new.append(each_trio.upper())

    score = 0
    for each_trio_upper in trio_list_new:
        if each_trio_upper in GAT_list:
            score += 4
        elif each_trio_upper in AAT_list:
            score += 3
        elif re.match(GAX, each_trio_upper):
            score += 2
        elif re.match(XAX, each_trio_upper):
            score += 1
        else:
            print('each_trio_upper error', each_trio_upper)

    return(score)


if __name__ == '__main__':

    file_trio = r'10_minmismatch_newspacer.csv'
    with open(file_trio, 'r') as f_trio:
        data_trio = pd.read_csv(f_trio)

        for i in data_trio.index:
            # print(i)
            nc = data_trio.at[i, 'Accession Number']

            trio_seq = data_trio.at[i, 'DnaA-trio']
            trio_seq_list = re.findall(r'.{3}', trio_seq)
            trio_score = get_trio_score(trio_seq_list)
            data_trio.at[i, 'trio score'] = trio_score

        data_trio.to_csv('11_score_trio.csv', index=False)
