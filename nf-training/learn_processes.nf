/*
Input values
The val qualifier allows you to receive data of any type as input. It can be accessed in the process script using the 
specified input name, as shown in the following example


num = Channel.of( 1, 2, 3)

process BASICEXAMPLE {
    debug true

    input:
    val x

    script:
    """
    echo process job $x
    """
}

workflow {
    myrun = BASICEXAMPLE(num)
}

Input files
The path qualifier allows the handling of file values in the process execution context

reads = Channel.fromPath( 'data/ggal/*.fq' )

process FOO {
    debug true

    input:
    path 'sample.fastq'

    script:

    """
    ls sample.fastq
    """
}

workflow{
    result = FOO(reads)
}

Input can also be defined using a variable. For example:


reads = Channel.fromPath( 'data/ggal/*.fq' )

process FOO {
    debug true

    input:
    path sample

    script:

    """
    ls $sample
    """
}

workflow{
    result = FOO(reads)
}

The same syntax is able to handle more than one input file in the same execution and only requires changing the channel composition


reads = Channel.fromPath( 'data/ggal/*.fq' )

process FOO {
    debug true

    input:
    path sample

    script:

    """
    ls -lh $sample
    """
}

workflow{
    result = FOO(reads.collect())
}
*/

/* Exercise: Write a script that concatenates all files ending with *_1.fastq into a new file and print out the first 20 lines


params.reads = 'data/ggal/*_1.fq' 
channel
    .fromPath( params.reads )
    .set{ read_ch }

process CONCATENATE {
    debug true
    tag "Concatenating files"


    input:
    path '*'

    output:
    path 'top_10'

    script:

    """
    cat * > total.csv
    head -n 20 total.csv > top_10
    """
}

workflow{
    concat_ch=CONCATENATE(read_ch.collect())
    concat_ch.view()
}

Write a process that is executed for each read file matching the pattern data/ggal/*_1.fq and use the same data/ggal/transcriptome.fa in each execution.


params.transcriptome_file = "$baseDir/data/ggal/transcriptome.fa"
params.reads = "data/ggal/*_1.fq"
channel
    .fromPath( params.reads )
    .set{ read_ch }


process COMMAND {
  tag "Run Command"

  input:
  path reads
  path transcriptome

  output:
  path result

  script:
   """
   echo your_command $reads $transcriptome > result
   """
}

workflow {
  concat_ch = COMMAND(read_ch, params.transcriptome_file)
  concat_ch.view()
}
*/

/* input repeaters
the each qualifier allows you to repeat the execution of a process for each item in a collection 
every time new data is received. for the mode, I can use e.g ['flye', 'raven', 'canu']


sequences = Channel.fromPath('data/prots/*.tfa')
methods = ['regular', 'espresso', 'psicoffee']

process ALIGNSEQUENCES {
  debug true

  input:
  path seq
  each mode

  script:
  """
  echo t_coffee -in $seq -mode $mode
  """
}

workflow {
  ALIGNSEQUENCES(sequences, methods)
}

Exercise: 
Extend the previous example so a task is executed for each read file matching the pattern data/ggal/*_1.fq and repeat the same task with both salmon and kallisto.


params.transcriptome_file = "$baseDir/data/ggal/transcriptome.fa"
params.reads = "$baseDir/data/ggal/*_1.fq"
methods = ['salmon', 'kallisto']
channel
    .fromPath( params.reads )
    .set{ read_ch }


process COMMAND {
  tag "Run Command"

  input:
  path reads
  path transcriptome
  each mode

  output:
  path result


  script:
  """
  echo $mode $reads $transcriptome > result
  """
}

workflow {
  concat_ch = COMMAND(read_ch, params.transcriptome_file, methods)
  concat_ch
    .view { "To run : ${it.text}" }
}

The val qualifier specifies a defined value in the script context. Values are frequently defined in the input and/or output declaration blocks, as shown in the following example:


methods = ['prot','dna', 'rna']

process FOO {

  input:
  val x

  output:
  val x

  script:
  """
  echo $x > file
  """
}

workflow {
  receiver_ch = FOO(Channel.of(methods))
  receiver_ch.view { "Received: $it" }
}

This printed out Received: [prot, dna, rna]

Chapter 6.3.2
The code below did not work


process RANDOMNUM {

    
    output:
    path 'result.txt'

    script:
    """
    echo $RANDOM > result.txt
    """
}


workflow {
  receiver_ch = RANDOMNUM()
  receiver_ch.view { "Received: " + it.text }
}

Multiple output files

When an output file name contains a wildcard character (* or ?) it is interpreted as a glob path matcher. This allows us to capture multiple files into a list object and output them as a sole emission.


process SPLITLETTERS {

    output:
    path 'chunk_*'

    """
    printf 'Hola' | split -b 1 - chunk_
    """
}

workflow {
    letters = SPLITLETTERS()
    letters
        .flatMap()
        .view { "File: ${it.name} => ${it.text}" }
}

without the .flatMap() the results look like this
File: [chunk_aa, chunk_ab, chunk_ac, chunk_ad] => [H, o, l, a]
with the .flatMap() the results look like this:
File: chunk_aa => H
File: chunk_ab => o
File: chunk_ac => l
File: chunk_ad => a

6.3.4: Dynamic output file names. I will need this for flye
When an output file name needs to be expressed dynamically, it is possible to define it using a dynamic string that references values defined in the input declaration block or in the script global context


species = ['cat','dog', 'sloth']
sequences = ['AGATAG','ATGCTCT', 'ATCCCAA']

Channel.fromList(species)
       .set { species_ch }

process ALIGN {

  input:
  val x
  val seq

  output:
  path "${x}.aln"

  script:
  """
  echo align -in $seq > ${x}.aln
  """
}

workflow {
  genomes = ALIGN( species_ch, sequences )
  genomes.view()
}

This produces e.g /mnt/c/Users/Afolayan/Documents/Courses/nextflow_training/nf-training-public/nf-training/work/60/966e6a178a6bb3d3c65616d6e7c1ef/dog.aln
within it, you find e.g., align -in [AGATAG, ATGCTCT, ATCCCAA]

6.3.5: Composite input and output


reads_ch = Channel.fromFilePairs('data/ggal/*_{1,2}.fq')

process FOO {

  input:
    tuple val(sample_id), path(sample_id)

  output:
    tuple val(sample_id), path("sample.bam")

  script:
  """
    echo your_command_here --reads $sample_id > sample_id.bam
  """
}

workflow {
  bam_ch = FOO(reads_ch)
  bam_ch.view()
}

the difference here is that the output remains a literal sample_id.bam for all files within the work dir
note that at the input section, the path is stated as path(sample_id)


reads_ch = Channel.fromFilePairs('data/ggal/*_{1,2}.fq')

process FOO {

  input:
    tuple val(sample_id), path(sample_files)

  output:
    tuple val(sample_id), path("${sample_id}.bam")

  script:
  """
    echo your_command_here --reads $sample_id > ${sample_id}.bam
  """
}

workflow {
  bam_ch = FOO(reads_ch)
  bam_ch.view()
}

the difference here is that the output remains a literal sample_id.bam for all files within the work dir
note that at the input section, the path is stated as path(sample_files) and the output script output name.bam file varies in name

6.4:  using the when delaration (functions like ifelse)
*/
/*
params.dbtype = 'nr'
params.prot = 'data/prots/*.tfa'
proteins = Channel.fromPath(params.prot)

process FIND {
    debug true

    input:
    path fasta
    val type

    when:
    fasta.name =~ /^BB11.* && type == 'nr' 

//    script:
//    """
//    echo blastp -query $fasta -db nr
//    """
//}

/*
workflow {
  result = FIND(proteins, params.dbtype)
}

6.5. Directives: allow the definition of optional settings that affect the execution of the current process without affecting the semantic of the task itself


process FOO {
  cpus 2
  memory 1.GB
  container 'image/name'

  script:
  """
  echo your_command --this --that
  """
}

6.6: Organize output, using the publishDir to publish the output


params.outdir = 'my-results'
params.prot = 'data/prots/*.tfa'
proteins = Channel.fromPath(params.prot)


process BLASTSEQ {
    publishDir "$params.outdir/bam_files", mode: 'copy'

    input:
    path fasta

    output:
    path ('*.txt')

    script:
    """
    echo blastp $fasta > ${fasta}_result.txt
    """
}

workflow {
  blast_ch = BLASTSEQ(proteins)
  blast_ch.view()
}
6.6.2: Manage semantic sub-directories


params.reads = 'data/reads/*_{1,2}.fq.gz'
params.outdir = 'my-results'

samples_ch = Channel.fromFilePairs(params.reads, flat: true)

process FOO {
  publishDir "$params.outdir/$sampleId/", pattern: '*.fq'
  publishDir "$params.outdir/$sampleId/counts", pattern: "*_counts.txt"
  publishDir "$params.outdir/$sampleId/outlooks", pattern: '*_outlook.txt'

  input:
    tuple val(sampleId), path('sample1.fq.gz'), path('sample2.fq.gz')

  output:
    path "*"

  script:
  """
    < sample1.fq.gz zcat > sample1.fq
    < sample2.fq.gz zcat > sample2.fq

    awk '{s++}END{print s/4}' sample1.fq > sample1_counts.txt
    awk '{s++}END{print s/4}' sample2.fq > sample2_counts.txt

    head -n 50 sample1.fq > sample1_outlook.txt
    head -n 50 sample2.fq > sample2_outlook.txt
  """
}

workflow {
  out_channel = FOO(samples_ch)
}

*/

/* Chapter 7.2.2
Use fromPath to create a channel emitting the fastq files matching the pattern data/ggal/*.fq, then use map to return a pair containing the file name and the path itself, and finally, use view to print the resulting channel
*/
/*
Channel
    .fromPath('data/ggal/*.fq')
    .map { it -> [it.name, it] }
    .view{ it, filepath -> "$it is found in $filepath"}

 e.g of answers, gut_1.fq is found in /mnt/c/Users/Afolayan/Documents/Courses/nextflow_training/nf-training-public/nf-training/data/ggal/gut_1.fq 
 */
/*
Channel
    .of(1,2,3)
    .collect()
    .view()

The result of the collect operator is a value channel

7.2.6 groupTuple
Exercise:
Use fromPath to create a channel emitting all of the files in the folder data/meta/, then use a map to associate the baseName prefix to each file. Finally, group all files that have the same common prefix

 
Channel
    .fromPath('data/meta/*')
    .map{ file -> tuple(file.baseName, file) }
    .groupTuple()
    .view { baseName, file -> "> $baseName : $file" }
    

*/

/*
8.9 for loop syntax
for (int i = 0; i <3; i++) {
   println("Hello World $i")
}

Iteration over list objects is also possible using the syntax below:
*/

list = ['a','b','c']

for( String elem : list ) {
  println elem
}


