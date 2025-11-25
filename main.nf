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
    path "functions.tsv", emit: func_table
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

process makeReport {
    publishDir params.out_dir, mode: "copy"
    label "r_env"

    input:
    path tax_table
    path func_table

    output:
    path "wf-faprotax-report.html"
    path "tables"

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(htmlwidgets)
    library(htmltools)
    library(reactable)

    get_table <- function(path) {
        reactable(
            fread(path, sep="\\t"),
            filterable=TRUE, compact=TRUE, striped=TRUE, bordered=TRUE,
            resizable=TRUE, pagination=FALSE, highlight=TRUE
        )
    }

    dir.create("tables", showWarnings = FALSE)

    tax_table <- get_table("${tax_table}")
    saveWidget(tax_table, "tables/tax_table.html", selfcontained = FALSE, libdir = "lib")
    func_table <- get_table("${func_table}")
    saveWidget(func_table, "tables/func_table.html", selfcontained = FALSE, libdir = "lib")
    page = tagList(
        h1("Report"),
        h2("Taxonomy"),
        tags\$iframe(
            src="tables/tax_table.html",
            frameBorder = "0",
            width="100%",
            height="600"
        ),
        h2("Functions"),
        tags\$iframe(
            src="tables/func_table.html",
            frameBorder = "0",
            width="100%",
            height="600"
        ),
    )
    save_html(page, "wf-faprotax-report.html")
    """
}

workflow {
    WorkflowMain.initialise(workflow, params, log)

    tax_table = file(params.table, checkIfExists: true)
    if (tax_table.isDirectory()) {
        ch_table = files("${tax_table}/*.tsv", checkIfExists: true)
        tax_table = joinTables(ch_table)
    }
    groups = file(params.groups, checkIfExists: true)
    convertToBiom(tax_table)
    faprotax(convertToBiom.out, groups)
    makeReport(tax_table, faprotax.out.func_table)
}
