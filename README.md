# clustering_commands

implementation of clustering algorithms with C++

# clustering_nonparametric_bayes_2d

This is the only command in this repository. This is an implementation of

https://github.com/Ma-sa-ue/practice/blob/master/machine%20learning(python)/nonparabayes/DPmixture(CRP).ipynb

with small modification. This command classifies two dimensional data on Re^2 with the assumption that they are based on a Dirichlet process Gaussian mixture model.

## input data format

There is an example file for input in "data" directory.

    $ cat data/example_input_2d
    1 1
    1 2
    2 1
    4 9
    4 8
    9 9
    10 9
    8 7
    9 10

The clustering\_nonparametric\_bayes\_2d command receives the data from its standard input and outputs it with cluster ids shown at the third column of the output. The numbers of the ids are scattered but they tells which records belong to each cluster.

    $ cat data/example_input_2d | ./clustering_nonparametric_bayes_2d 2> /dev/null
    1 1 2
    1 2 2
    2 1 2
    4 9 0
    4 8 0
    9 9 3
    10 9 3
    8 7 3
    9 10 3

The standard error tells the progression of clustering.

    $ cat data/example_input_2d | ./clustering_nonparametric_bayes_2d > /dev/null
    sweep 48
    0.33,0.83 cov: 0.02,0.02 num: 2
    0.037,0.037 cov: 0.02,0.02 num: 3
    0.89,0.86 cov: 0.02,0.02 num: 4
    ----
    sweep 49
    0.33,0.83 cov: 0.02,0.02 num: 2
    0.037,0.037 cov: 0.02,0.02 num: 3
    0.89,0.86 cov: 0.02,0.02 num: 4
    ----
    
