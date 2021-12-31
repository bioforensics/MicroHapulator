# USC panel

> <small>Example configuration files for a panel published by [de la Puente, Phillips, el al. (FSI Genetics, 2019)](https://doi.org/10.1016/j.fsigen.2019.102213).</small>

Use the Python interpreter to get marker IDs.

```python
>>> import microhapdb
>>> import pandas as pd
>>> m = microhapdb.markers
>>> usc_panel = m[(m.Source == "10.1016/j.fsigen.2019.102213") & ~(pd.isna(m.Ae))]
>>> usc_panel[["Name"]].to_csv("usc-panel.txt", index=False, header=False)
```

Then use shell commands for the rest.

```bash
$ # Grab reference sequences for the panel
$ microhapdb marker --format=fasta --delta=25 --min-length=200 --panel=usc-panel.txt > usc-refr.fasta
$
$ # Grab marker definitions, i.e., the SNP offsets for each marker
$ microhapdb marker --format=offsets --delta=25 --min-length=200 --panel=usc-panel.txt > usc-defn.tsv
$
$ # Grab haplotype frequency estimates for "Iberian Population in Spain" calculated from
$ # 1000 Genomes Project data.
$ microhapdb frequency --format=mhpl8r --panel=usc-panel.txt --population=IBS > usc-freq-esp.tsv
```
