#!/usr/bin/env python3
from collections import defaultdict
import pandas
import sys
param_table = pandas.read_csv(sys.stdin, sep='\t')
params = defaultdict(list)
for n, values in param_table.iterrows():
    sampleid = values['ID']
    params[sampleid].append(values)
for sampleid, param_lists in sorted(params.items()):
    outfile = sampleid + '.fastq.gz'
    if len(param_lists) == 1:
        gtfile = str(sampleid) + '.genotype.bed'
        cmd = 'mhpl8r sim --out {of:s} --genotype {gt:s} --panel usa --hapseed {hs:d} --num-reads 250000 --seq-seed {ss:d} --set-threads 1 {pop1:s} {pop2:s}'.format(
            of=outfile, gt=gtfile, hs=param_lists[0]['HaploSeed'], ss=param_lists[0]['SeqSeed'],
            pop1=param_lists[0]['MaternalHaploPop'], pop2=param_lists[0]['PaternalHaploPop']
        )
    else:
        indivstr = ' '.join(['--indiv {} {}'.format(p['MaternalHaploPop'], p['PaternalHaploPop']) for p in param_lists])
        propstr = ' '.join([str(p['Proportion']) for p in param_lists])
        haplostr = ' '.join([str(p['HaploSeed']) for p in param_lists])
        seqstr = ' '.join([str(p['SeqSeed']) for p in param_lists])
        cmd = 'mhpl8r mixture {indiv:s} --panel usa --out {of:s} --proporions {prop:s} --hap-seeds {hs:s} --seq-seeds {ss:s} --seq-threads 1'.format(
           indiv=indivstr, of=outfile, prop=propstr, hs=haplostr, ss=seqstr
        )
    print(cmd)
