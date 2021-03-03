import pandas as pd
import regex as rx
import re
import itertools
import copy


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
    regex_box_str_1 = '(' + regex_box_str + ')'
    return(regex_box_str_1)


def box1_box2(result, motif):
    result = result.reset_index(drop=True)
    for i in result.index:
        box1 = result.at[i, 'box1 in trio']
        if str(result.at[i, 'GC_rich_region']) == 'nan':
            gc_region = ''
        else:
            gc_region = result.at[i, 'GC_rich_region']
        seq_query = box1[-1] + gc_region
        for mis in range(5):
            boxRegex = re.compile(box_regex(motif, mis), re.I)
            if re.search(boxRegex, seq_query):
                m = re.search(boxRegex, seq_query)
                result.at[i, 'box2 in trio'] = m.group()
                result.at[i, 'box2 mismatch'] = mis
                result.at[i, 'seq between box1 and box2'] = seq_query[1:m.start()]
                result.at[i, 'GC_rich_region'] = seq_query[m.end():]
                result.at[i, 'the length of GC_rich_region'] = len(
                    seq_query[m.end():])

                if m.start() == 0:
                    result.at[i, 'box1 in trio'] = box1[:-1]

                break
            else:
                continue
    return(result)


def box1_box2_specific_box(result):
    result = result.reset_index(drop=True)
    for i in result.index:

        motif = result.at[i, 'box motif']
        box1 = result.at[i, 'box1 in trio']
        if str(result.at[i, 'GC_rich_region']) == 'nan':
            gc_region = ''
        else:
            gc_region = result.at[i, 'GC_rich_region']
        seq_query = box1[-1] + gc_region
        for mis in range(5):
            boxRegex = re.compile(box_regex(motif, mis), re.I)
            if re.search(boxRegex, seq_query):
                m = re.search(boxRegex, seq_query)
                result.at[i, 'box2 in trio'] = m.group()
                result.at[i, 'box2 mismatch'] = mis
                result.at[i, 'seq between box1 and box2'] = seq_query[1:m.start()]
                result.at[i, 'GC_rich_region'] = seq_query[m.end():]
                result.at[i, 'the length of GC_rich_region'] = len(
                    seq_query[m.end():])

                if m.start() == 0:
                    result.at[i, 'box1 in trio'] = box1[:-1]

                break
            else:
                continue
    return(result)


if __name__ == '__main__':

    file = r'4_concat_delete_repeat.csv'
    with open(file, 'r') as f:
        data = pd.read_csv(f)

    groups = data.groupby(by=['box motif'])

    data_standard = pd.DataFrame()
    data_specific = pd.DataFrame()
    for box, group in groups:
        if box == 'AATAGGTGT':
            data_standard = pd.concat([group, data_standard])
        else:
            data_specific = pd.concat([group, data_specific])

    data_standard = data_standard.reset_index(drop=True)
    data_specific = data_specific.reset_index(drop=True)

    result_standard = box1_box2(data_standard, 'AATAGGTGT')
    result_specific = box1_box2_specific_box(data_specific)

    result = pd.concat([result_standard, result_specific])

    result.to_csv('5_box1_box2.csv', index=False)
