process foo {
  debug true
  '''
  env | egrep 'ALPHA|BETA'
  '''
}


workflow {
    foo()
}
