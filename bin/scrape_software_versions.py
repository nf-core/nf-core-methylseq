#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'nf-core/methylseq': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'Bismark genomePrep': ['v_bismark_genome_preparation.txt', r"Bismark Genome Preparation Version: v(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'Cutadapt': ['v_cutadapt.txt', r"(\S+)"],
    'Trim Galore!': ['v_trim_galore.txt', r"version (\S+)"],
    'Bismark': ['v_bismark.txt', r"Bismark Version: v(\S+)"],
    'Bismark Deduplication': ['v_deduplicate_bismark.txt', r"Deduplicator Version: v(\S+)"],
    'Bismark methXtract': ['v_bismark_methylation_extractor.txt', r"Bismark Extractor Version: v(\S+)"],
    'Bismark Report': ['v_bismark2report.txt', r"bismark2report version: v(\S+)"],
    'Bismark Summary': ['v_bismark2summary.txt', r"bismark2summary version: (\S+)"],
    'HISAT2': ['v_hisat2.txt', r"version (\S+)"],
    'Samtools': ['v_samtools.txt', r"samtools (\S+)"],
    'BWA': ['v_bwa.txt', r"Version: (\S+)"],
    'bwa-meth': ['v_bwameth.txt', r"bwa-meth\.py (\S+)"],
    'Picard MarkDuplicates': ['v_picard_markdups.txt', r"([\d\.]+)"],
    'Picard CreateSequenceDictionary': ['v_picard_createseqdict.txt', r"([\d\.]+)"],
    'Picard CollectInsertSizeMetrics': ['v_picard_collectinssize.txt', r"([\d\.]+)"],
    'Picard CollectGcBiasMetrics': ['v_picard_collectgcbias.txt', r"([\d\.]+)"],
    'samblaster': ['v_samblaster.txt', r"samblaster: Version (\S+)"],
    'biscuit': ['v_biscuit.txt', r"Version: (\S+)"],
    'bcftools': ['v_bcftools.txt', r"bcftools (\S+)"],
    'bedtools': ['v_bedtools.txt', r"bedtools (\S+)"],
    'parallel': ['v_parallel.txt', r"GNU parallel (\S+)"],
    'gawk': ['v_gawk.txt', r"GNU Awk (\S+)"],
    'MethylDackel': ['v_methyldackel.txt', r"(.+)"],
    'Qualimap': ['v_qualimap.txt', r"QualiMap v.(\S+)"],
    'Preseq': ['v_preseq.txt', r"Version: (\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
    
}
results = OrderedDict()
results['nf-core/methylseq'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['Bismark genomePrep'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['Cutadapt'] = '<span style="color:#999999;\">N/A</span>'
results['Trim Galore!'] = '<span style="color:#999999;\">N/A</span>'
results['Bismark'] = '<span style="color:#999999;\">N/A</span>'
results['Bismark Deduplication'] = '<span style="color:#999999;\">N/A</span>'
results['Bismark methXtract'] = '<span style="color:#999999;\">N/A</span>'
results['Bismark Report'] = '<span style="color:#999999;\">N/A</span>'
results['Bismark Summary'] = '<span style="color:#999999;\">N/A</span>'
results['Samtools'] = '<span style="color:#999999;\">N/A</span>'
results['BWA'] = '<span style="color:#999999;\">N/A</span>'
results['bwa-meth'] = '<span style="color:#999999;\">N/A</span>'
results['Picard MarkDuplicates'] = '<span style="color:#999999;\">N/A</span>'
results['MethylDackel'] = '<span style="color:#999999;\">N/A</span>'
results['Qualimap'] = '<span style="color:#999999;\">N/A</span>'
results['Preseq'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'
results['samblaster'] = '<span style="color:#999999;\">N/A</span>'
results['biscuit'] = '<span style="color:#999999;\">N/A</span>'
results['fastasort'] = '<span style="color:#999999;\">N/A</span>'
results['Picard CreateSequenceDictionary'] = '<span style="color:#999999;\">N/A</span>'
results['Picard CollectInsertSizeMetrics'] = '<span style="color:#999999;\">N/A</span>'
results['Picard CollectGcBiasMetrics'] = '<span style="color:#999999;\">N/A</span>'
results['bcftools'] = '<span style="color:#999999;\">N/A</span>'
results['bedtools'] = '<span style="color:#999999;\">N/A</span>'
results['gawk'] = '<span style="color:#999999;\">N/A</span>'
results['parallel'] = '<span style="color:#999999;\">N/A</span>'


# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
    except IOError:
        results[k] = False

# Remove software set to false in results
for k in list(results):
    if not results[k]:
        del results[k]

# Dump to YAML

print(
    """
id: 'software_versions'
section_name: 'nf-core/methylseq Software Versions'
section_href: 'https://github.com/nf-core/methylseq'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
"""
)
for k, v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k, v))
print("    </dl>")

# Write out regexes as csv file:
with open("software_versions.csv", "w") as f:
    for k, v in results.items():
        f.write("{}\t{}\n".format(k, v))
