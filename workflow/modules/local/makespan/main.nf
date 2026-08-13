process PREPARE_GENE_SPANS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.11'
        : 'quay.io/biocontainers/python:3.11'}"

    input:
    path ref_gtf

    output:
    path "gene_spans.gtf", emit: gtf
    tuple val("${task.process}"), val('python'), eval("python3 --version"), topic: versions, emit: versions_preparegenespans

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3
    import re, sys, gzip

    gene_id_re = re.compile(r'gene_id "([^"]+)"')
    gene_name_re = re.compile(r'gene_name "([^"]+)"')

    def opener(path):
        return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

    genes = {}
    in_path, out_path = "${ref_gtf}", "gene_spans.gtf"

    with opener(in_path) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\\n').split('\\t')
            if len(parts) < 9:
                continue
            chrom, _, _, start_s, end_s, _, strand, _, attrs = parts
            try:
                start, end = int(start_s), int(end_s)
            except ValueError:
                continue
            gid_m = gene_id_re.search(attrs)
            if not gid_m:
                continue
            gene_id = gid_m.group(1)
            gname_m = gene_name_re.search(attrs)
            gene_name = gname_m.group(1) if gname_m else gene_id
            if gene_id not in genes:
                genes[gene_id] = {'chrom': chrom, 'start': start, 'end': end,
                                   'strand': strand, 'name': gene_name, 'id': gene_id}
            else:
                g = genes[gene_id]
                if chrom != g['chrom']:
                    continue
                g['start'] = min(g['start'], start)
                g['end'] = max(g['end'], end)

    records = sorted(genes.values(), key=lambda g: (g['chrom'], g['start']))

    with open(out_path, "w") as out:
        for g in records:
            out.write(f"{g['chrom']}\\tsynthetic\\tgene\\t{g['start']}\\t{g['end']}\\t.\\t{g['strand']}\\t.\\t"
                       f'gene_id "{g["id"]}"; gene_name "{g["name"]}";\\n')

    print(f"genes written: {len(records)}", file=sys.stderr)
    """

    stub:
    """
    touch gene_spans.gtf
    """
}
