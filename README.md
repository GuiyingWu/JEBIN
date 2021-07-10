# User Guide of JEBIN

## Introduction

This is the JEBIN algorithm developed for learning the consensus representations of genes by joint embedding of multiple bipitite networks. 


## Install
Before running JEBIN on Linux system, packages [GSL] (https://www.gnu.org/software/gsl/)  and [Eigen] (http://eigen.tuxfamily.org/index.php?title=Main_Page) are required to be installed. Users must modify the package paths before using the makefile to compile the code.


## Usage
```
./JEBIN_linux -consensus_nodes_file consensus_nodes_list.txt -network_filenames_file network_filenames_list.txt -num_network 2 -output_directory output/ -binary 0 -size 200 -negative 5 -samples 200 -threads 10 -gamma 1 -rho 0.025
```

- -consensus_nodes_file: 
                the consensus gene set
- -network_filenames_file: 
                the file contains all the bipartite network filenames (using absolute paths)
- -num_network: 
                the number of networks
- -output_directory: 
                the output directory of all kinds of the nodes' vectors
- -binary: 
                save the resulting vectors in binary mode; default is 0 (off)
- -size: 
                the dimension of the embedding vectors; default is 200
- -negative: 
                the number of negative examples used in negative sampling; default is 5
- -samples: 
                the total number of training samples (*Million)
- -threads: 
                the total number of threads used; the default is 1
- -gamma: 
                the regularizing coefficient; the default is 1.0
- -rho: 
                the starting value of the learning rate; the default is 0.025




## Input 

#### consensus_nodes_file 

This file contains the union of genes of all the datasets to be integrated (Gene Entrez ID is recommended). An example is shown bellow: 
```
100
1000
10001
10006
10007
10008
10009
100093630
```

#### network_filenames_file

This file contains the absolute paths of all the bipartite network files with each constructed from one gene expression dataset. An example is shown below:
```
/home/gywu/multiset/data/bulkHCC/edgelist_Dataset1_HCCDB1.txt
/home/gywu/multiset/data/bulkHCC/edgelist_Dataset2_HCCDB13.txt
/home/gywu/multiset/data/bulkHCC/edgelist_Dataset3_HCCDB15.txt
/home/gywu/multiset/data/bulkHCC/edgelist_Dataset4_HCCDB17.txt
/home/gywu/multiset/data/bulkHCC/edgelist_Dataset5_HCCDB18.txt
/home/gywu/multiset/data/bulkHCC/edgelist_Dataset6_HCCDB3.txt
/home/gywu/multiset/data/bulkHCC/edgelist_Dataset7_HCCDB4.txt
/home/gywu/multiset/data/bulkHCC/edgelist_Dataset8_HCCDB6.txt
```

Each bipartite network file contains the edges between genes (first column) and samples (second column), the last column is the normalized gene expression value. An example is shown below:
```
100	HCCDB-1.S3	7.5599
100	HCCDB-1.S5	8.3189
100	HCCDB-1.S6	7.4908
100	HCCDB-1.S8	7.1945
100	HCCDB-1.S10	8.104
100	HCCDB-1.S12	7.9868
100	HCCDB-1.S14	8.5352
100	HCCDB-1.S16	7.6303
```


The "/data" folder contains three examples of the input data. 

"/scHCC" folder contains the single-cell RNA-seq data of HCC (gene filtered), which is in "rds" format.

"/output" folder contains three examples of the output results of JEBIN:
- "output_u_consensus.txt" contains the consensus representation vectors for genes across all networks.
- "output_u_net1.txt" contains the dataset-specific representation vectors for genes in the first input network.
- "output_v_net1.txt" contains the dataset-specific representation vectors for samples in the first input network.



## Contact
```
Guiying Wu (email: wgy14@mails.tsinghua.edu.cn)
```


## Citation
```
@article{wu2020new,
title={New gene association measures by joint network embedding of multiple gene expression datasets},
author={Wu, Guiying and Li, Xiangyu and Guo, Wenbo and Wei, Zheng and Hu, Tao and Gu, Jin},
journal={bioRxiv},
year={2020}}
```
