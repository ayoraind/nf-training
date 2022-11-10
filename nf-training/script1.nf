params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.multiqc = "$projectDir/multiqc"
params.outdir = "nf-training/results"

/*
* println "reads: $params.reads"
* println "reads: $params.reads"
*/

log.info """\
         R N A S E Q - N F   P I P E L I N E    
         ===================================
         transcriptome: ${params.transcriptome_file}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

/*
 * define the `index` process that creates a binary index
 * given the transcriptome file
 */
process INDEX {
  input:
  path transcriptome

  output:
  path 'salmon_index'

  script:
  """
  salmon index --threads $task.cpus -t $transcriptome -i salmon_index
  """
}

workflow {
    index_ch = INDEX(params.transcriptome_file)
    index_ch.view{ it }
}

/* the view method throws in the following filepath
/mnt/c/Users/Afolayan/Documents/Courses/nextflow_training/nf-training-public/nf-training/work/34/696395d817ceed2b96b840d71f1e44/salmon_index
*/