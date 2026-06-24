# CovGain: An efficient and cost-effective filtering methodology to improve reads recovery in CpG-context large-scale DNA methylation studies

## Background

Population-level epigenomics provides crucial insights. However, state of the art techniques like Enzymatic Methyl Sequencing (EM-seq) remain financially and computationally demanding. Thus, the optimization of existing methods is of great interest for the field. MethylKit, a widely used analysis software, applies a strand-specific minimum coverage threshold and then, very frequently, merges both strands CpG information to increase coverage. This logic inadvertently discards CpG sites that fail the threshold individually but would pass if merged beforehand. It also reduces final coverage for loci passing due to only one strand, leading to severe information loss.

## Results

Here, we present CovGain, an R extension for methylKit that implements a combined-strand coverage filtering strategy. Validated using two European seabass (Dicentrarchus labrax) EM-seq datasets representing different tissue types (liver vs. caudal fin) and library strategies (pooling vs. individual samples) across a sequencing depth gradient (2×, 5×, 10×, 15× and 30×), CovGain significantly outperformed the standard methylKit methodology in CpG recovery. Notably, at 15×, CovGain recovered ~77% of genomic CpG sites per sample, a yield the standard method only achieves at ~25×, representing a 22.1% increase. Furthermore, while the standard method at 15× virtually failed to construct a CpG matrix for population analysis, CovGain successfully consolidated ~4-5 million common CpG sites without altering global methylation patterns and capturing ~60% of the Differentially Methylated Cytosines (DMCs) identified by the 30× standard method.

## Conclusion

CovGain allows to obtain valuable DNAm data at 15× sequencing depth reducing per-sample costs by approximately 30% and halving bioinformatic storage and processing demands. CovGain is a cost-effective and accessible software that allows to prioritize larger sample sizes maximizing statistical power in epigenomics, which is of particular interest for population-scale studies

