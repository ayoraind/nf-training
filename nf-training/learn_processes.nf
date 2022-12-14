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

*/

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