#!/usr/bin/env bash
# workflow/contract/test/download_test_data.sh
# Downloads HyperTRIBE validation dataset (GSE102814)
# Drosophila dm6 reference genome + 4 SRA samples

set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


mkdir -p "${DIR}/ref_genome" "${DIR}/raw_reads"

echo "Downloading dm6 reference genome..."
[ -f "${DIR}/ref_genome/dm6.fa" ] || {
    wget -P "${DIR}/ref_genome" https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
    gunzip "${DIR}/ref_genome/dm6.fa.gz"
}

echo "Decompressing GTF annotation..."
[ -f "${DIR}/ref_genome/genes.gtf" ] || {
    gunzip -k "${DIR}/ref_genome/genes.gtf.gz"
}

echo "Downloading SRA samples..."
for acc in SRR5944748 SRR5944749 SRR6426146 SRR5944750; do
    [ -f "${DIR}/raw_reads/${acc}.fastq" ] || {
        prefetch "${acc}" -O "${DIR}/raw_reads"
        fastq-dump "${DIR}/raw_reads/${acc}/${acc}.sra" -O "${DIR}/raw_reads"
    }
done

echo "Done. Runnning the test profile with:"
echo "nextflow run workflow/main.nf -profile test,docker"
