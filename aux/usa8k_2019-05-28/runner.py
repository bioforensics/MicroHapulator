#!/usr/bin/env python3
from collections import defaultdict
import pandas
import sys


def pop_str(pop1, pop2):
    if pop1 == pop2:
        return pop1
    return pop1 + ' ' + pop2


param_table = pandas.read_csv(sys.stdin, sep='\t')
params = defaultdict(list)
for n, values in param_table.iterrows():
    sampleid = values['ID']
    params[sampleid].append(values)
for sampleid, param_lists in sorted(params.items()):
    outfile = sampleid + '.fastq.gz'
    if len(param_lists) == 1:
        gtfile = str(sampleid) + '.genotype.bed'
        pops = pop_str(param_lists[0]['MaternalHaploPop'], param_lists[0]['PaternalHaploPop'])
        cmd = 'mhpl8r sim --out {of:s} --genotype {gt:s} --panel usa --hap-seed {hs:d} --num-reads 250000 --seq-seed {ss:d} --seq-threads 1 {pop:s}'.format(
            of=outfile, gt=gtfile, hs=param_lists[0]['HaploSeed'], ss=param_lists[0]['SeqSeed'],
            pop=pops
        )
    else:
        indivstr = ' '.join(['--indiv ' + pop_str(p['MaternalHaploPop'], p['PaternalHaploPop']) for p in param_lists])
        propstr = ' '.join([str(p['Proportion']) for p in param_lists])
        haplostr = ' '.join([str(p['HaploSeed']) for p in param_lists])
        seqstr = ' '.join([str(p['SeqSeed']) for p in param_lists])
        cmd = 'mhpl8r mixture {indiv:s} --panel usa --out {of:s} --genotype {sid:s}.bed --proportions {prop:s} --hap-seeds {hs:s} --seq-seeds {ss:s} --seq-threads 1'.format(
           indiv=indivstr, of=outfile, prop=propstr, hs=haplostr, ss=seqstr, sid=sampleid
        )
    print(cmd)
