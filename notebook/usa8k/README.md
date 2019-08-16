# USA 8k

```bash
# Produce parameters for mock genotypes and simulated sample sequencing
snakemake prep

# Create mock genotypes, simulate sample sequencing, map reads, infer genotypes,
#   and perform sample matching tests
snakemake --config seqthreads=8 --cores $THREADS all_sample_matching_tests
```
