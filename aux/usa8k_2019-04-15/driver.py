#!/usr/bin/env python3
import gen_sim_params

simulators = [
    gen_sim_params.gen_params(4000, 'usa8k-simple-', ncontrib=1, seed=4294195554),
    gen_sim_params.gen_params(400, 'usa8k-mixture-2-even-', ncontrib=2, even=True, seed=2762204106),
    gen_sim_params.gen_params(400, 'usa8k-mixture-3-even-', ncontrib=3, even=True, seed=2124137260),
    gen_sim_params.gen_params(400, 'usa8k-mixture-4-even-', ncontrib=4, even=True, seed=778059273),
    gen_sim_params.gen_params(400, 'usa8k-mixture-5-even-', ncontrib=5, even=True, seed=2972194569),
    gen_sim_params.gen_params(400, 'usa8k-mixture-6-even-', ncontrib=6, even=True, seed=2570983910),
    gen_sim_params.gen_params(400, 'usa8k-mixture-2-uneven-', ncontrib=2, seed=481612327),
    gen_sim_params.gen_params(400, 'usa8k-mixture-3-uneven-', ncontrib=3, seed=2019626428),
    gen_sim_params.gen_params(400, 'usa8k-mixture-4-uneven-', ncontrib=4, seed=3191594573),
    gen_sim_params.gen_params(400, 'usa8k-mixture-5-uneven-', ncontrib=5, seed=597241784),
    gen_sim_params.gen_params(400, 'usa8k-mixture-6-uneven-', ncontrib=6, seed=2822645581),
]

columns = ['ID', 'Label', 'MaternalHaploPop', 'PaternalHaploPop', 'HaploSeed', 'SeqSeed', 'Proportion']
print(*columns, sep='\t')
for simulator in simulators:
    for data in simulator:
        print(*data, sep='\t')
