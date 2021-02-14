# METAL-TPP

Metal Extraction Triggered Agitation Logged by Thermal Proteome Profiling (METAL-TPP) is designed for global discovery of MBPs in proteomes, which operates by extracting metals from MBPs with chelators and logging the resulting structural perturbation of MBPs with thermal proteome profiling. Protocol and scripts for analyzing  the TMT quantification data is provided in this repo.

Contact: chuwang@pku.edu.cn, xinzeng@pku.edu.cn, wendao@pku.edu.cn



##  Data pre-processing

1. Prepare TMT data file (intensities of the six reporter ions for each protein from Maxquant),  for example, PBS_proteingroup.txt and EDTA_proteingroup.txt for each replicate.

2. Calculate correction factors for each temperature

```bash
iter_fit_rawdata.py PBS_proteingroup.txt EDTA_proteingroup.txt
```

3. For each replicate, two sets of six correction factors were obtained.



## Fitting Melting curves

1. For each TMT data file, six correction factors are needed

2. Recalibrate the TMT data

```bash
parse_TMT.py [PBS/EDTA]_proteingroup.txt f1 f2 f3 f4 f5 f6 > [PBS/EDTA]_calibrated.txt
```

