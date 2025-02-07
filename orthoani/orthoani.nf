workflow {
    // define the parameters
    input = Channel
        .fromPath("./input/*.LargeContigs.fna")
        .map{tuple( it.name.split('.LargeContigs.fna')[0], it )}
    // run the process
    orthoani(input)
}

process orthoani {
  module 'Java/17.0.6'
  publishDir "output", mode: "copy"
  input:
    tuple val(accession), path(contigs)
  output:
    tuple val(accession), path("${accession}.log")
  errorStrategy 'ignore'
  script:
  """
  java -jar OAT_cmd.jar -num_threads ${task.cpus} -fasta1 "./ref/ref.fna" -fasta2 "./input/${accession}.LargeContigs.fna" -blastplus_dir ./blast_executable/ncbi-blast/bin/ >> "${accession}.log"
  """
}
