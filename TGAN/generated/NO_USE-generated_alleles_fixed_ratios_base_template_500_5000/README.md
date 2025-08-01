## Dataset Description
* ILS peak heights: 3000
* Base template: 500, 5000

```
degradation_settings <- list(
  list(shape = 2.5, scale = 1e-3),
  list(shape = 3.5, scale = 2e-3)
)
```

```
allele_freqs_file <- system.file("extdata","FBI_extended_Cauc_022024.csv", package = "simDNAmixtures")
```