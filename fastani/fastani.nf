workflow {
    // define the parameters
    input = Channel
        .fromPath("./input/*.LargeContigs.fna")
        .map{tuple( it.name.split('.LargeContigs.fna')[0], it )}
    // run the process
    fastani(input)
}

process fastani {
  container "fastani_1.34--h4dfc31f_0.sif"
  publishDir "output", mode: "copy"
  input:
    tuple val(accession), path(contigs)
  output:
    tuple val(accession), path("${accession}.log")
  errorStrategy 'ignore'
  script:
  """
  fastANI -t ${task.cpus} -q ./ref/ref.fna -r "${contigs}" -o "${accession}.log"
  """
}
