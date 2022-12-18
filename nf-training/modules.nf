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
    stdout emit: upper

    """
    cat $y | tr '[a-z]' '[A-Z]' 
    """
}

// define one or more parameters and functions

params.foo = 'Hello'
params.bar = 'world!'

def SAYHELLO() {
    println "$params.foo $params.bar"
}
