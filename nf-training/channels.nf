/*ch = Channel.of(1,2,3)
ch.view()


Channel
  .of(1..23, 'X', 'Y')
  .view()

/* channel from list


list = ['hello', 'world']

Channel
  .fromList(list)
  .view()

/* channel from path
*/
/*
Channel
    .fromPath( './data/meta/*.csv' )
    .view()


/* channel from another path


Channel
    .fromPath( './data/ggal/*.fq' , hidden:true)
    .view()
*/
/*output of the code below is a tuple with two elements, which are: sample name, and list of matching fastqs e.g.
[lung, [/mnt/c/Users/Afolayan/Documents/Courses/nextflow_training/nf-training-public/nf-training/data/ggal/lung_1.fq, /mnt/c/Users/Afolayan/Documents/Courses/nextflow_training/nf-training-public/nf-training/data/ggal/lung_2.fq]]


Channel
  .fromFilePairs('./data/ggal/*_{1,2}.fq')
  .view()


/*output of the code below is a tuple with three elements, which are: sample name, and each individual element (not grouped as a list as in the previous code). E.g.
[lung, /mnt/c/Users/Afolayan/Documents/Courses/nextflow_training/nf-training-public/nf-training/data/ggal/lung_1.fq, /mnt/c/Users/Afolayan/Documents/Courses/nextflow_training/nf-training-public/nf-training/data/ggal/lung_2.fq]

Channel
  .fromFilePairs('./data/ggal/*_{1,2}.fq', flat:true)
  .view()
*/

/*fromSRA. This code prints contents of an NCBI project ID


params.ncbi_api_key = '92d6a5cd2a5bafed87caafbba1b6cce7f10a'

Channel
  .fromSRA(['SRP073307'], apiKey: params.ncbi_api_key)
  .view()

/*fromSRA. This code prints multiple accession ID

params.ncbi_api_key = '92d6a5cd2a5bafed87caafbba1b6cce7f10a'

ids = ['ERR908507', 'ERR908506', 'ERR908505']
Channel
  .fromSRA(ids, apiKey: params.ncbi_api_key)
  .view()

/* now run fastqc
process FASTQC {
    container = 'biocontainers/fastqc:v0.11.5'
    tag "FASTQC on $sample_id"
...
*/

/*run fastqc on SRA accessions
*/

params.ncbi_api_key = '92d6a5cd2a5bafed87caafbba1b6cce7f10a'

params.accession = ['ERR908507', 'ERR908506']

process fastqc {
    container = 'biocontainers/fastqc:v0.11.5'
    tag "FASTQC on $sample_id"
    input:
    tuple val(sample_id), path(reads_file)

    output:
    path("fastqc_${sample_id}_logs")

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads_file}
    """
}

workflow {
    reads = Channel.fromSRA(params.accession, apiKey: params.ncbi_api_key)
    fastqc(reads)
}