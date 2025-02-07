workflow {
    // define the parameters
    input = Channel
        .fromPath("./input/*.LargeContigs.fna")
        .map{tuple( it.name.split('.LargeContigs.fna')[0], it )}
    // run the process
    skani(input)
}

process skani {
  container "skani_0.2.0--h4ac6f70_0.sif"
  publishDir "output", mode: "copy"
  input:
    tuple val(accession), path(contigs)
  output:
    tuple val(accession), path("${accession}.log")
  errorStrategy 'ignore'
  script:
  """
  skani dist -t ${task.cpus} -q ./ref/ref.fna -r "${contigs}" -o "${accession}.log"
  """
}
