#!/usr/bin/env bash
# contract/test/download_test_data.sh
# Downloads HyperTRIBE validation dataset (GSE102814)
# Drosophila dm6 reference + 4 SRA samples

set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Downloading dm6 reference genome..."
wget -P "${DIR}/ref_genome" https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
wget -P "${DIR}/ref_genome" https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/genes/dm6.ensGene.gtf.gz
gunzip "${DIR}/ref_genome/dm6.fa.gz"
gunzip "${DIR}/ref_genome/dm6.ensGene.gtf.gz"

echo "Downloading SRA samples..."
for acc in SRR5944748 SRR5944749 SRR6426146 SRR5944750; do
    prefetch "${acc}" -O "${DIR}/raw_reads"
    fastq-dump "${DIR}/raw_reads/${acc}/${acc}.sra" -O "${DIR}/raw_reads"
done

echo "Done. Run the test profile with:"
echo "nextflow run workflow/main.nf -profile test,docker"
