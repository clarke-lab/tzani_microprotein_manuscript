 
# CHO Cell Microproteins 

[![DOI](https://zenodo.org/badge/449655379.svg)](https://zenodo.org/badge/latestdoi/449655379)

The code contained in this repositority enable the reproduction of the results of:

Castro-Rivadeneyra *et. al* 2023. **Annotation of the non-canonical translatome reveals that CHO cell microproteins are a new class of mAb drug product impurity**

The publication is freely availiable here: xxxxxxx&nbsp;

**Abstract:**
<p style='text-align: justify;'>
Mass spectrometry (MS) has emerged as a powerful approach for the detection of Chinese hamster ovary (CHO) cell protein impurities in antibody drug products. The incomplete annotation of the Chinese hamster genome, however, limits the coverage of MS-based host cell protein (HCP) analysis.</p> &nbsp;

<p style='text-align: justify;'>
Chinese hamster ovary (CHO) cells are used to produce almost 90% of therapeutic monoclonal antibodies (mAbs). The annotation of non-canonical translation events in these cellular factories remains incomplete, limiting not only our ability to study CHO cell biology but also detect host cell protein (HCP) contaminants in the final mAb drug product. We utilised ribosome footprint profiling (Ribo-seq) to identify novel open reading frames (ORFs) including N-terminal extensions and thousands of short ORFs (sORFs) predicted to encode microproteins. Mass spectrometry-based HCP analysis of four commercial mAb drug products using the extended protein sequence database revealed the presence of microprotein impurities for the first time. We also show that microprotein abundance varies with growth phase and can be affected by the cell culture environment. In addition, our work provides a vital resource to facilitate future studies of non-canonical translation as well as the regulation of protein synthesis in CHO cell lines.
</p>
&nbsp;



# Preparation


## Mamba Environment

### 1. Linux software
```bash
mamba env create -f microprotein_process_env.yaml --name microprotein_process_env
```

### 2. R and R packages
```bash
mamba env create -f microprotein_r_env.yaml --name microprotein_r_env
```

## Conda Environment
### 1. MetaMorpheus 
```
https://github.com/smith-chem-wisc/MetaMorpheus/wiki/Getting-Started#test-conda-installation-linux-macos-windows
```

## Docker
### 1. ORF-RATER
```bash
docker pull clarkelab/orfrater
```

## Analysis
## 1. Process the NGS data
```bash
./scripts/process_ngs_data.sh
```

## 2. Novel ORF identification
```bash
./scripts/identify_novel_orfs.sh
```

## 3. Translational Efficency Analysis

```bash
./scripts/run_differential_translation.sh
```

## 4. Mass Spectrometry Data analysis
```bash
./scripts/run_proteomics.sh
```

## Manuscript
```bash
./scripts/prepare_manuscript.sh
```