// By default, the process command is interpreted as a Bash script. However, any other scripting language can be used by simply starting the script with the corresponding Shebang declaration. For example:

// process PYSTUFF {

//   script:
//   """
//   #!/usr/bin/env python

//   x = 'Hello'
//   y = 'world!'
//   print ("%s - %s" % (x,y))
//   """
// }

// workflow {
//   PYSTUFF()
// }

// Script parameters (params) can be defined dynamically using variable values. For example:

// params.data = 'World'

// process FOO {

//   script:
//   """
//   echo Hello $params.data
//   """
// }

// workflow {
//   FOO()
// }

// The process script can also be defined in a completely dynamic manner using an if statement or any other expression for evaluating a string value. For example:
// params.compress = 'gzip'
// params.file2compress = "$baseDir/data/ggal/transcriptome.fa"

// process FOO {

//   input:
//   path file

//   script:
//   if( params.compress == 'gzip' )
//     """
//     gzip -c $file > ${file}.gz
//     """
//   else if( params.compress == 'bzip2' )
//     """
//     bzip2 -c $file > ${file}.bz2
//     """
//   else
//     throw new IllegalArgumentException("Unknown aligner $params.compress")
// }

// workflow {
//   FOO(params.file2compress)
// }

// The val qualifier allows you to receive data of any type as input.

// num = Channel.of( 1, 2, 3 )

// process BASICEXAMPLE {
//   debug true

//   input:
//   val x

//   script:
//   """
//   echo process job $x
//   """
// }

// workflow {
//   myrun = BASICEXAMPLE(num)
// }

// The path qualifier allows the handling of file values in the process execution context. This means that Nextflow will stage it in the process execution directory, and it can be accessed in the script by using the name specified in the input declaration.

// reads = Channel.fromPath( 'data/ggal/*.fq' )

// process FOO {
//   debug true

//   input:
//   path 'sample.fastq'

//   script:
//   """
//   ls sample.fastq
//   """
// }

// workflow {
//   result = FOO(reads)
// }

// The input file name can also be defined using a variable reference as shown below:

// reads = Channel.fromPath( 'data/ggal/*.fq' )

// process FOO {
//   debug true

//   input:
//   path sample

//   script:
//   """
//   ls  $sample
//   """
// }

// workflow {
//   result = FOO(reads)
// }

// The same syntax is also able to handle more than one input file in the same execution and only requires changing the channel composition.

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

workflow {
  FOO(reads.collect())
}

// Write a script that creates a channel containing all read files matching the pattern data/ggal/*_1.fq followed by a process that concatenates them into a single file and prints the first 20 lines.