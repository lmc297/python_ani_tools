workflow {
    // define the parameters
    input = Channel
        .fromPath("./input/*.LargeContigs.fna")
        .map{tuple( it.name.split('.LargeContigs.fna')[0], it )}
    // run the process
    pyfastani(input)
}

process pyfastani {
  container "pyfastani.sif"
  publishDir "output", mode: "copy"
  input:
    tuple val(accession), path(contigs)
  output:
    tuple val(accession), path("${accession}.log")
  errorStrategy 'ignore'
  script:
  """
  python ./run_pyfastani.py ./ref/ref.fna "${contigs}" ${task.cpus} "${accession}.log"
  """
}
