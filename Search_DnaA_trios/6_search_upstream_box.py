import pandas as pd
import re
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import copy
import itertools


def box_list(motif):
    motif_list = []
    if '[' in motif or ']' in motif:
        index_list = []
        for i in range(len(motif)):
            if motif[i] == '[' or motif[i] == ']':
                index_list.append(i)

        if index_list[0] != 0:
            index = index_list[0]
            motif_list.append(motif[:index])

        for j in range(len(index_list)-1):
            if int(j) % 2 == 0:
                index1 = index_list[j]
                index2 = index_list[j+1]+1
                motif_list.append(motif[index1:index2])
            elif int(j) % 2 == 1:
                index1 = index_list[j]+1
                index2 = index_list[j+1]
                motif_list.append(motif[index1:index2])
            else:
                print('box_list error', j)

        if index_list[-1] < (len(motif)-1):
            index = index_list[-1]+1
            motif_list.append(motif[index:])

        if '' in motif_list:
            motif_list.remove('')

        motif_list_result = []
        for n in motif_list:
            if '[' in n or ']' in n:
                motif_list_result.append(n)
            else:
                motif_list_result.extend(list(n))
    else:
        motif_list_result = list(motif)
    return(motif_list_result)


def box_regex(motif, mismatch):
    regex_set = set()
    tp_constant = box_list(motif)

    items = itertools.combinations(range(len(tp_constant)), mismatch)

    for item in items:
        tp = copy.deepcopy(tp_constant)
        for i in range(mismatch):
            tp[item[i]] = '.'
        regex_set.add(''.join(tp))
    regex_box_str = r'|'.join(regex_set)

    return(regex_box_str)


def GetSeqLen(now_spacer_len):
    SeqLen = 16-now_spacer_len
    if SeqLen < 0:
        return(0)
    else:
        return(SeqLen)


if __name__ == '__main__':

    file_trio = r'5_box1_box2.csv'
    with open(file_trio, 'r') as f1:
        data_trio = pd.read_csv(f1)

    file_doric = r'DoriC_oric_sequences.csv'
    with open(file_doric, 'r') as f2:
        data_doric = pd.read_csv(f2)
        data_doric.index = data_doric['DoriC AC']

    for i in data_trio.index:

        if str(data_trio.at[i, 'box2 in trio']) == 'nan':

            ori = data_trio.at[i, 'ORI']
            trio = data_trio.at[i, 'box_gc_trio_seq']

            if str(data_trio.at[i, 'GC_rich_region']) == 'nan':
                now_spacer_len = len(data_trio.at[i, 'box1 in trio'])
            else:
                now_spacer_len = len(
                    data_trio.at[i, 'box1 in trio'] + data_trio.at[i, 'GC_rich_region'])

            motif = data_trio.at[i, 'box motif']
            motif_list = box_list(motif)
            motif_len = len(motif_list)

            oric_seq = data_doric.at[ori, 'oric seq'].upper()

            oric_seq = Seq(oric_seq, IUPAC.unambiguous_dna)
            seq1 = str(oric_seq)[::-1]
            strand1 = 1
            seq2 = str(oric_seq.complement())
            strand2 = 2

            pattern = re.compile(trio, re.I)

            if re.search(pattern, seq1):
                match = re.search(pattern, seq1)

                box0_search_seq_start = match.start() - GetSeqLen(now_spacer_len) - len(motif_list)
                box0_search_seq = seq1[box0_search_seq_start:match.start(
                )] + trio[0]

                motif = data_trio.at[i, 'box motif']

                for mismatch in range(5):
                    box_regex_str = box_regex(motif, mismatch)
                    box_regex_expression = re.compile(box_regex_str)
                    if re.search(box_regex_expression, box0_search_seq):
                        match_box = re.search(
                            box_regex_expression, box0_search_seq)
                        if match_box.end() == len(box0_search_seq):
                            print(i, ori, match_box.end(),
                                  len(box0_search_seq))
                            data_trio.at[i, 'box_0'] = match_box.group()[:-1]
                            data_trio.at[i, 'box0_seq_box1'] = box0_search_seq[match_box.end(
                            ):-1]
                            data_trio.at[i, 'box0 mismatch'] = mismatch
                            break
                        else:
                            data_trio.at[i, 'box_0'] = match_box.group()
                            data_trio.at[i, 'box0_seq_box1'] = box0_search_seq[match_box.end(
                            ):-1]
                            data_trio.at[i, 'box0 mismatch'] = mismatch
                            break

            elif re.search(pattern, seq2):
                match = re.search(pattern, seq2)

                box0_search_seq_start = match.start() - GetSeqLen(now_spacer_len) - len(motif_list)
                box0_search_seq = seq2[box0_search_seq_start:match.start(
                )] + trio[0]

                motif = data_trio.at[i, 'box motif']

                for mismatch in range(5):
                    box_regex_str = box_regex(motif, mismatch)
                    box_regex_expression = re.compile(box_regex_str)
                    if re.search(box_regex_expression, box0_search_seq):
                        match_box = re.search(
                            box_regex_expression, box0_search_seq)
                        if match_box.end() == len(box0_search_seq):
                            print(i, ori, match_box.end(),
                                  len(box0_search_seq))
                            data_trio.at[i, 'box_0'] = match_box.group()[:-1]
                            data_trio.at[i, 'box0_seq_box1'] = box0_search_seq[match_box.end(
                            ):-1]
                            data_trio.at[i, 'box0 mismatch'] = mismatch
                            break
                        else:
                            data_trio.at[i, 'box_0'] = match_box.group()
                            data_trio.at[i, 'box0_seq_box1'] = box0_search_seq[match_box.end(
                            ):-1]
                            data_trio.at[i, 'box0 mismatch'] = mismatch
                            break
            else:
                print(ori)

    data_trio.to_csv('6_search_upstream_box.csv', index=False)
