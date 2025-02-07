workflow {
    // define the parameters
    input = Channel
        .fromPath("./input/*.LargeContigs.fna")
        .map{tuple( it.name.split('.LargeContigs.fna')[0], it )}
    // run the process
    pyorthoani(input)
}

process pyorthoani {
  container "./pyorthoani.sif"
  publishDir "output", mode: "copy"
  input:
    tuple val(accession), path(contigs)
  output:
    tuple val(accession), path("${accession}.log")
  errorStrategy 'ignore'
  script:
  """
  orthoani -q "./ref/ref.fna" -r "${contigs}" -j ${task.cpus} >> "${accession}.log"
  """
}
