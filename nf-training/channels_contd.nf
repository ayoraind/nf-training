/*The splitText operator allows you to split multi-line strings or text file items, emitted by a source channel into chunks containing n lines, which will be emitted by the resulting channel


Channel
  .fromPath('data/meta/random.txt') 
  .splitText()                      
  .view()

/*You can define the number of lines in each chunk by using the parameter by, as shown in the following example:
The subscribe operator permits execution of user defined functions each time a new value is emitted by the source channel.


Channel
  .fromPath('data/meta/random.txt')
  .splitText( by: 2 )
  .subscribe {
    print it;
    print "--- end of the chunk ---\n"
  }

/*An optional closure can be specified in order to transform the text chunks produced by the operator. The following example shows how to split text files into chunks of 10 lines and transform them into capital letters:


Channel
  .fromPath('data/meta/random.txt')
  .splitText( by: 10 ) { it.toUpperCase() }
  .view()

/*You can also make counts for each line:

count=0

Channel
  .fromPath('data/meta/random.txt')
  .splitText()
  .view { "${count++}: ${it.toUpperCase().trim()}" }


/*Finally, you can also use the operator on plain files (outside of the channel context):

def f = file('data/meta/random.txt')
def lines = f.splitText()
def count=0
for( String row : lines ) {
  log.info "${count++} ${row.toUpperCase()}"
}

/*apply the splitCsv operator to a channel emitting a CSV formatted text files or text entries. For example, to view only the first and fourth columns


Channel
  .fromPath("data/meta/patients_1.csv")
  .splitCsv()
  // row is a list object
  .view { row -> "${row[0]},${row[3]}" }

/*When the CSV begins with a header line defining the column names, you can specify the parameter header: true which allows you to reference each value by its column name, as shown in the following example:


Channel
  .fromPath("data/meta/patients_1.csv")
  .splitCsv(header: true)
  // row is a list object
  .view { row -> "${row.patient_id},${row.num_samples}" }

/*Alternatively, you can provide custom header names by specifying a list of strings in the header parameter as shown below:


Channel
  .fromPath("data/meta/patients_1.csv")
  .splitCsv(header: ['col1', 'col2', 'col3', 'col4', 'col5'] )
  // row is a list object
  .view { row -> "${row.col1},${row.col4}" }


//You can also process multiple csv files at the same time:
*/

// Channel
//   .fromPath("data/meta/patients_*.csv") // <-- just use a pattern
//   .splitCsv(header:true)
//   .view { row -> "${row.patient_id}\t${row.num_samples}" }


// //Finally, you can also operate on csv files outside the channel context:
// def f = file('data/meta/patients_1.csv')
//   def lines = f.splitCsv()
//   for( List row : lines ) {
//     log.info "${row[0]} -- ${row[2]}"
//   }

//Task: Try inputting fastq reads into the RNA-seq workflow from earlier using .splitCSV

//step1: Add a csv text file containing the following, as an example input with the name "fastq.csv":

        //gut,/workspace/nf-training-public/nf-training/data/ggal/gut_1.fq,/workspace/nf-training-public/nf-training/data/ggal/gut_2.fq
        //Then replace the input channel for the reads in script7.nf. Changing the workflow line to a splitCSV channel factory:

// workflow {
//     Channel
//         .fromPath("fastq.csv") // <-- just use a pattern
//         .splitCsv()
//         .view () { row -> "${row[0]},${row[1]},${row[2]}" }
//         .set { read_pairs_ch }

//     index_ch = INDEX(params.transcriptome_file)
//     quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)
//     fastqc_ch = FASTQC(read_pairs_ch)
//     MULTIQC(quant_ch.mix(fastqc_ch).collect())
// }

//step2: Finally, change the cardinality of the processes that use the input data. For example, for the quantification process, to
/* process QUANTIFICATION {
  tag "$sample_id"

  input:
  path salmon_index
  tuple val(sample_id), path(reads1), path(reads2)

  output:
  path sample_id, emit: quant_ch

  script:
  """
  salmon quant --threads $task.cpus --libType=U -i $salmon_index -1 ${reads1} -2 ${reads2} -o $sample_id
  """
}
*/

// do same for the FASTQ process

/*process FASTQC {
  tag "FASTQC on $sample_id"

  input:
  tuple val(sample_id), path(reads1), path(reads2)

  output:
  path "fastqc_${sample_id}_logs"

  script:
  """
  mkdir fastqc_${sample_id}_logs
  fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads1} ${reads2}
  """
}
*/

// Parsing tsv files works in a similar way, simply add the sep:'\t' option in the splitCsv context:

// Channel
//   .fromPath("data/meta/regions.tsv", checkIfExists:true)
//   // use `sep` option to parse TAB separated files
//   .splitCsv(sep:'\t')
//   // row is a list object
//   .view()

// skip header

// Channel
//   .fromPath("data/meta/regions.tsv", checkIfExists:true)
//   // use `sep` option to parse TAB separated files
//   .splitCsv(sep:'\t', header:true)
//   // row is a list object
//   .view { row -> "${row.patient_id}" }


// parsing json files

// import groovy.json.JsonSlurper

// def f = file('data/meta/regions.json')
// def records = new JsonSlurper().parse(f)


// for( def entry : records ) {
//   log.info "$entry.patient_id -- $entry.feature"
// }

// parsing yaml files

import org.yaml.snakeyaml.Yaml

def f = file('data/meta/regions.yml')
def records = new Yaml().load(f)


for( def entry : records ) {
  log.info "$entry.patient_id -- $entry.feature"
}