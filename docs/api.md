# Python API Reference

**NOTE**: The MicroHapulator Python API is under [Semantic Versioning](https://semver.org/).
In brief, this means that every stable version of the MicroHapulator software is assigned a version number, and that any changes to the software's behavior or interface require the software version number to be updated in prescribed and predictable ways.


## Haplotype calling

```{eval-rst}
.. autofunction:: microhapulator.api.type
```

```{eval-rst}
.. autofunction:: microhapulator.profile.TypingResult.filter
```


## Analysis, QA/QC, and interpretation

```{eval-rst}
.. autofunction:: microhapulator.api.read_length_dist
```

```{eval-rst}
.. autofunction:: microhapulator.api.interlocus_balance
```

```{eval-rst}
.. autofunction:: microhapulator.api.plot_haplotype_calls
```

```{eval-rst}
.. autofunction:: microhapulator.api.heterozygote_balance
```

```{eval-rst}
.. autofunction:: microhapulator.api.contrib
```

```{eval-rst}
.. autofunction:: microhapulator.api.prob
```

```{eval-rst}
.. autofunction:: microhapulator.api.diff
```

```{eval-rst}
.. autofunction:: microhapulator.api.dist
```

```{eval-rst}
.. autofunction:: microhapulator.api.contain
```

## Simulation

```{eval-rst}
.. autofunction:: microhapulator.api.sim
```

```{eval-rst}
.. autofunction:: microhapulator.profile.SimulatedProfile.merge
```

```{eval-rst}
.. autofunction:: microhapulator.profile.Profile.unite
```

```{eval-rst}
.. autofunction:: microhapulator.api.seq
```
