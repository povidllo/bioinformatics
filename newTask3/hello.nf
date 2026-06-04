process hello {
    output:
    stdout

    script:
    """
    echo "Hello Nextflow"
    """
}

workflow {
    hello | view
}
