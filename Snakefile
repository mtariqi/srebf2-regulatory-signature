# ==============================================================================
# Project     : Chromatin-state context of SREBF2 knockout-responsive genes in HeLa cells: an exploratory integrative epigenomic analysis
# File        : Snakefile
# Purpose     : SREBF2 Regulatory Signature
# Author      : Md Tariqul Islam
# Institution : Northeastern University | Bioinformatics Department
# License     : MIT (see LICENSE)
# ==============================================================================

# SREBF2 Regulatory Signature
# Phase 12 reproducibility workflow
#
# Uses the existing Phase 12E regulatory-signal matrix.
# Statistical scripts are not modified.

SIGNAL = "results/phase12/SREBF2_HELA_REGULATORY_SIGNAL.tsv"

rule all:
    input:
        "reports/phase12/phase12_signal_qc.log",
        "results/phase12/PHASE12_EXPLORATORY_COMPARISONS.tsv",
        "results/phase12/PHASE12_ADJUSTED_MODELS.tsv",
        "results/phase12/PHASE12_DIRECTION_MODELS.tsv",
        "results/phase12/PHASE12_ROBUSTNESS_MODELS.tsv",
        "results/phase12/PHASE12_VIF.tsv",
        "results/phase12/PHASE12_EXPRESSION_STRATIFIED_MODELS.tsv",
        "results/phase12/PHASE12_EXPRESSION_STRATA_COUNTS.tsv"

rule phase12_signal_qc:
    input:
        SIGNAL,
        "scripts/phase12/phase12_signal_qc.py"
    output:
        "reports/phase12/phase12_signal_qc.log"
    shell:
        """
        python3 scripts/phase12/phase12_signal_qc.py \
          > {output} 2>&1
        """

rule phase12_exploratory:
    input:
        SIGNAL,
        "scripts/phase12/phase12_exploratory_comparison.py"
    output:
        "results/phase12/PHASE12_EXPLORATORY_COMPARISONS.tsv"
    log:
        "reports/phase12/phase12_exploratory_comparison.log"
    shell:
        """
        python3 scripts/phase12/phase12_exploratory_comparison.py \
          > {log} 2>&1
        """

rule phase12_adjusted:
    input:
        SIGNAL,
        "scripts/phase12/phase12_adjusted_models.py"
    output:
        "results/phase12/PHASE12_ADJUSTED_MODELS.tsv"
    log:
        "reports/phase12/phase12_adjusted_models.log"
    shell:
        """
        python3 scripts/phase12/phase12_adjusted_models.py \
          > {log} 2>&1
        """

rule phase12_direction:
    input:
        SIGNAL,
        "scripts/phase12/phase12_direction_models.py"
    output:
        "results/phase12/PHASE12_DIRECTION_MODELS.tsv"
    log:
        "reports/phase12/phase12_direction_models.log"
    shell:
        """
        python3 scripts/phase12/phase12_direction_models.py \
          > {log} 2>&1
        """

rule phase12_robustness:
    input:
        SIGNAL,
        "scripts/phase12/phase12_robustness.py"
    output:
        "results/phase12/PHASE12_ROBUSTNESS_MODELS.tsv",
        "results/phase12/PHASE12_VIF.tsv"
    log:
        "reports/phase12/phase12_robustness.log"
    shell:
        """
        python3 scripts/phase12/phase12_robustness.py \
          > {log} 2>&1
        """

rule phase12_expression_stratified:
    input:
        SIGNAL,
        "scripts/phase12/phase12_expression_stratified.py"
    output:
        "results/phase12/PHASE12_EXPRESSION_STRATIFIED_MODELS.tsv",
        "results/phase12/PHASE12_EXPRESSION_STRATA_COUNTS.tsv"
    log:
        "reports/phase12/phase12_expression_stratified.log"
    shell:
        """
        python3 scripts/phase12/phase12_expression_stratified.py \
          > {log} 2>&1
        """
