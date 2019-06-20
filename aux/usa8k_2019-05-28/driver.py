#!/usr/bin/env python3
from collections import defaultdict
import gen_sim_params
import numpy

simulators = [
    gen_sim_params.gen_params(4000, ncontrib=1, seed=4294195554),
    gen_sim_params.gen_params(400, ncontrib=2, even=True, seed=2762204106),
    gen_sim_params.gen_params(400, ncontrib=3, even=True, seed=2124137260),
    gen_sim_params.gen_params(400, ncontrib=4, even=True, seed=778059273),
    gen_sim_params.gen_params(400, ncontrib=5, even=True, seed=2972194569),
    gen_sim_params.gen_params(400, ncontrib=6, even=True, seed=2570983910),
    gen_sim_params.gen_params(400, ncontrib=2, even=False, seed=481612327),
    gen_sim_params.gen_params(400, ncontrib=3, even=False, seed=2019626428),
    gen_sim_params.gen_params(400, ncontrib=4, even=False, seed=3191594573),
    gen_sim_params.gen_params(400, ncontrib=5, even=False, seed=597241784),
    gen_sim_params.gen_params(400, ncontrib=6, even=False, seed=2822645581),
]

samples = list()
samples_by_population = defaultdict(list)

with open('sample-contrib.tsv', 'w') as cout, open('sample-seq.tsv', 'w') as sout:
    columns = ['ID', 'Contributors', 'SeqSeeds', 'Proportions']
    print(*columns, sep='\t', file=sout)
    columns = ['ID', 'MaternalHaploPop', 'PaternalHaploPop', 'HaploSeed']
    print(*columns, sep='\t', file=cout)
    for simulator in simulators:
        for sample in simulator:
            if sample.ncontrib == 1:
                samples.append(sample)
                samples_by_population[sample.popstring].append(sample)
            print(sample.seqparams, file=sout)
            print(sample.simparams, file=cout)

    numpy.random.seed(3709894119)
    with open('sample-unite.tsv', 'w') as uout:
        columns = ['ID', 'MaternalID', 'PaternalID', 'Seed']
        print(*columns, sep='\t', file=uout)
        for _ in range(100):
            probe = numpy.random.choice(samples)
            target_pop = samples_by_population[probe.popstring]
            mom, dad = numpy.random.choice(target_pop, 2, replace=False)
            kidlabel = gen_sim_params.get_id()
            print(
                kidlabel, mom.contriblabel, dad.contriblabel,
                numpy.random.randint(2**32 - 1), sep='\t', file=uout
            )
            print(
                gen_sim_params.get_id(), kidlabel, numpy.random.randint(2**32 - 1),
                '1.0', file=sout
            )
