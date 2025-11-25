process joinTables {
    publishDir params.out_dir, mode: "copy"
    label "miller"

    input:
    path tables

    output:
    path "table.tsv"

    script:
    """
    mlr --tsv put '\$[[1]] = "taxonomy"' \\
        then unsparsify --fill-with 0 \\
        then stats1 -a sum -g taxonomy --fx '^taxonomy\$' \\
        then rename -r '_sum\$,' ${tables} > table.tsv
    """
}

process convertToBiom {
    publishDir params.out_dir, mode: "copy"
    label "faprotax"

    input:
    path table

    output:
    path "table.biom"

    script:
    """
    biom convert -i ${table} -o table.biom --to-hdf5
    """
}

process faprotax {
    publishDir params.out_dir, mode: "copy"
    label "faprotax"

    input:
    path biom
    path groups

    output:
    path "functions.tsv", emit: funct_table
    path "report.txt", emit: report

    script:
    options = []
    if (params.omit_unrepresented_groups) {
        options << "--omit_unrepresented_groups"
    }
    options = options.join(" ")
    groups = groups.name == "NO_FILE" ? "/app/FAPROTAX.txt" : groups
    """
    collapse_table.py \\
        -i ${biom} \\
        -g ${groups} \\
        ${options} \\
        -o functions.tsv \\
        -r report.txt
    """
}

workflow {
    WorkflowMain.initialise(workflow, params, log)

    table = file(params.table, checkIfExists: true)
    if (table.isDirectory()) {
        ch_table = files("${table}/*.tsv", checkIfExists: true)
        table = joinTables(ch_table)
    }
    groups = file(params.groups, checkIfExists: true)
    convertToBiom(table)
    faprotax(convertToBiom.out, groups)
}
