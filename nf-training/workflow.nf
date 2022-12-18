#!/usr/bin/env nextflow

params.greeting = 'Hello world!'

/*
include { SPLITLETTERS } from './modules.nf'
include { CONVERTTOUPPER } from './modules.nf'

// A workflow component can declare one or more input channels using the take statement
// When the take statement is used, the workflow definition needs to be declared within the main block.
workflow my_pipeline {
    take:
    greeting

    main:
    SPLITLETTERS(greeting)
    CONVERTTOUPPER(SPLITLETTERS.out.flatten())
    CONVERTTOUPPER.out.upper.view{ it }
}

workflow {
    my_pipeline(Channel.of(params.greeting))
}

*/ 

// 9.3.3: Calling named workflows
/*
Within a main.nf script (called hello.nf in our example) we also can have multiple workflows. In which case we may want to call a specific workflow when running the code. For this we use the entrypoint call -entry <workflow_name>
*/

include { SPLITLETTERS as SPLITLETTERS_one } from './modules.nf'
include { SPLITLETTERS as SPLITLETTERS_two } from './modules.nf'

include { CONVERTTOUPPER as CONVERTTOUPPER_one } from './modules.nf'
include { CONVERTTOUPPER as CONVERTTOUPPER_two } from './modules.nf'


workflow my_pipeline_one {
    letters_ch1 = SPLITLETTERS_one(params.greeting)
    results_ch1 = CONVERTTOUPPER_one(letters_ch1.flatten())
    results_ch1.view{ it }
}

workflow my_pipeline_two {
    letters_ch2 = SPLITLETTERS_two(params.greeting)
    results_ch2 = CONVERTTOUPPER_two(letters_ch2.flatten())
    results_ch2.view{ it }
}

workflow {
    my_pipeline_one(Channel.of(params.greeting))
    my_pipeline_two(Channel.of(params.greeting))
}