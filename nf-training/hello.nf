#!/usr/bin/env nextflow

// params.greeting = 'Hello world!'
// greeting_ch = Channel.of(params.greeting)

/*
process SPLITLETTERS {
    input:
    val x

    output:
    path 'chunk_*'

    """
    printf '$x' | split -b 6 - chunk_
    """
}

process CONVERTTOUPPER {
    input:
    path y

    output:
    stdout

    """
    cat $y | tr '[a-z]' '[A-Z]' 
    """
}

workflow {
    letters_ch = SPLITLETTERS(greeting_ch)
    results_ch = CONVERTTOUPPER(letters_ch.flatten())
    results_ch.view{ it }
}
*/

/*
include { SPLITLETTERS   } from './modules.nf'
include { CONVERTTOUPPER } from './modules.nf'
*/

/*
workflow {
    letters_ch = SPLITLETTERS(greeting_ch)
    results_ch = CONVERTTOUPPER(letters_ch.flatten())
    results_ch.view{ it }
}
*/
// seems like when you use emit, you have to use .out function in your workflow channel outputs
/*
workflow {
    greeting_ch = Channel.of(params.greeting)
    SPLITLETTERS(greeting_ch)
    CONVERTTOUPPER(SPLITLETTERS.out.flatten())
    CONVERTTOUPPER.out.upper.view{ it }
}
*/

/*
workflow {
    Channel.of(params.greeting) | SPLITLETTERS | flatten() | CONVERTTOUPPER | view
}
*/

/*
workflow my_pipeline {
    greeting_ch = Channel.of(params.greeting)
    SPLITLETTERS(greeting_ch)
    CONVERTTOUPPER(SPLITLETTERS.out.flatten())
    CONVERTTOUPPER.out.upper.view{ it }
}

workflow {
    my_pipeline()
}
*/

/*
params.foo = 'Hola'
params.bar = 'mundo!'

include { SAYHELLO } from './modules.nf'

workflow {
    SAYHELLO()
}

*/

// The addParams option can be used to extend the module parameters without affecting the external scope
params.foo = 'Hola'
params.bar = 'mundo!'

include { SAYHELLO } from './modules.nf' addParams(foo: 'Olá')

workflow {
    SAYHELLO()
}