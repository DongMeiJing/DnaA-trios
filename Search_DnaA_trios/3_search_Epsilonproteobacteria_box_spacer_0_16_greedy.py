import pandas as pd
import re
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import copy
import itertools
import regex as rx


def which_box(NC):
    with open(r'Epsilonproteobacteria_DnaA_Box_motif.csv', 'r') as f1:
        data_box = pd.read_csv(f1)
    nc_box_list = list(data_box['NC'])
    if NC in nc_box_list:
        index_box = data_box[data_box['NC'] == NC].index
        if len(index_box) == 1:
            motif = data_box.at[index_box[0], 'complementary box']
        else:
            print('which_box error:more than one box')
    else:
        motif = 'AATAGGTGT'
    return(motif)


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
    tp_constant = box_list(motif)[::-1]

    items = itertools.combinations(range(len(tp_constant)), mismatch)

    for item in items:
        tp = copy.deepcopy(tp_constant)
        for i in range(mismatch):
            tp[item[i]] = '.'
        regex_set.add(''.join(tp))
    regex_box_str = r'|'.join(regex_set)
    regex_box_str_1 = '(' + regex_box_str + ')'
    return(regex_box_str_1)


def GC_region_regex(max_gc_len):
    gc_regex = '([ACGT]{{0,{}}})'.format(max_gc_len)
    return(gc_regex)


def which_regex(motif, mis, max_gc_len):
    regex = rx.compile(
        '((.A.){3,})' + GC_region_regex(max_gc_len) + box_regex(motif, mis), re.I)
    return(regex)


def start_end(start_old, end_old, start_new, end_new, strand, oric_seq, maxPosition):
    if strand == 1:
        start = start_old + start_new
        end = start_old + end_new - 1
        if start_old > end_old:
            if start > maxPosition and end > maxPosition:
                start = start - maxPosition
                end = end - maxPosition
            if end > maxPosition and start <= maxPosition:
                end = end - maxPosition
        return(start, end)

    if strand == 2:
        end = end_old - end_new+1
        start = end_old - start_new
        if start_old > end_old:
            if start < 0 and end < 0:
                end = maxPosition-end_new+end_old+1
                start = maxPosition-start_new+end_old
            if end < 0 and start >= 0:
                end = maxPosition-end_new+end_old+1
        return(start, end)


def search_trio(pattern, seq, start_old, end_old, strand, NC, mis, motif, maxPosition, ORI):
    result = pd.DataFrame()
    if rx.search(pattern, seq):
        line = 0
        for everyseq in rx.finditer(pattern, seq, overlapped=True):
            find_seq = everyseq.group()

            start_new = everyseq.start()
            end_new = everyseq.end()
            start = start_end(start_old, end_old, start_new,
                              end_new, strand, seq, maxPosition)[0]
            end = start_end(start_old, end_old, start_new,
                            end_new, strand, seq, maxPosition)[1]

            result.at[line, 'NC'] = NC
            result.at[line, 'ORI'] = ORI
            result.at[line, 'box_gc_trio_seq'] = find_seq[::-1]
            result.at[line, 'the start of trio'] = start
            result.at[line, 'the end of trio'] = end
            result.at[line, 'box motif'] = motif
            result.at[line, 'box1 mismatch'] = mis
            result.at[line, 'strand'] = strand
            result.at[line, 'box1 in trio'] = everyseq.group(4)[::-1]
            result.at[line, 'GC_rich_region'] = everyseq.group(3)[::-1]
            result.at[line, 'the length of GC_rich_region'] = len(
                everyseq.group(3))
            result.at[line, 'trio'] = everyseq.group(1)[::-1]
            result.at[line, 'the number of trio'] = len(everyseq.group(1))/3
            result.at[line, 'chromosome size'] = maxPosition
            line += 1
        return(result)
    return(result)


def all_steps(seq, start, end, NC, motif, maxPosition, mis_max, max_gc_len, ORI):
    result = pd.DataFrame()

    seq = Seq(seq, IUPAC.unambiguous_dna)
    seq1 = str(seq)
    strand1 = 1
    seq2 = str(seq.reverse_complement())
    strand2 = 2

    for mis in range(mis_max+1):
        pattern = which_regex(motif, mis, max_gc_len)
        result1 = search_trio(pattern, seq1, start, end,
                              strand1, NC, mis, motif, maxPosition, ORI)
        result2 = search_trio(pattern, seq2, start, end,
                              strand2, NC, mis, motif, maxPosition, ORI)
        result = pd.concat([result, result1, result2], axis=0)
    return(result)


if __name__ == '__main__':

    result = pd.DataFrame()

    file_Epsilonproteobacteria = r'Epsilonproteobacteria_DnaA_Box_motif.csv'
    with open(file_Epsilonproteobacteria, 'r') as f_Epsilonproteobacteria:
        data_Epsilonproteobacteria = pd.read_csv(f_Epsilonproteobacteria)
        nc_Epsilonproteobacteria = list(data_Epsilonproteobacteria['NC'])

    file_doric = r'DoriC_oric_sequences.csv'
    with open(file_doric, 'r') as f_doric:
        data_doric = pd.read_csv(f_doric)

        for i in data_doric.index:
            NC = data_doric.at[i, 'NC']

            if NC in nc_Epsilonproteobacteria:
                print(NC)
                ori = data_doric.at[i, 'DoriC AC']
                oric_seq = data_doric.at[i, 'oric seq'].upper()
                start = int(data_doric.at[i, 'the start of oric'])
                end = int(data_doric.at[i, 'the end of oric'])
                maxPosition = int(data_doric.at[i, 'Genome Length'])

                max_gc_len = 16
                motif = which_box(NC)
                motif_list = box_list(motif)
                if len(motif_list) == 9:
                    mis_max = 2
                elif len(motif_list) > 9:
                    print(motif)
                    mis_max = 3
                else:
                    print(motif)

                each_result = all_steps(
                    oric_seq, start, end, NC, motif, maxPosition, mis_max, max_gc_len, ori)
                result = pd.concat([result, each_result])

        result.to_csv(
            '3_search_Epsilonproteobacteria_box_spacer_0_16_greedy.csv', index=False)
