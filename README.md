# User Guide of JEBIN

## Introduction

This is the JEBIN algorithm developed for learning the consensus representations of genes by joint embedding of multiple bipitite networks. 


## Install
Before running JEBIN on Linux system, packages [GSL] (https://www.gnu.org/software/gsl/)  and [Eigen] (http://eigen.tuxfamily.org/index.php?title=Main_Page) are required to be installed.


## Usage
```
./JEBIN_linux -consensus_nodes_file consensus_nodes_list.txt -network_filenames_file network_filenames_list.txt -num_network 2 -output_directory output/ -binary 0 -size 200 -negative 5 -samples 200 -threads 10 -gamma 1 -rho 0.025
```

- -consensus_nodes_file <file>
                The consensus nodes set
- -network_filenames_file <file>
                All the bipartite network files' names (absolute paths)
- -num_network <int>
                the number of networks
- -output_directory <file>
                The output directory of all kinds of the nodes' vectors
- -binary <int>
                Save the resulting vectors in binary moded; default is 0 (off)
- -size <int>
                Set size of word vectors; default is 100
- -negative <int>
                Number of negative examples; default is 5, common values are 5 - 10 (0 = not used)
- -samples <float>
                Set the number of training samples as <float>Million; default is 1 million
- -threads <int>
                Use <int> threads (default 1)
- -gamma <float>
                Set the regularizing coefficient; default is 1.0
- -rho <float>
                Set the starting learning rate; default is 0.025




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
/home/gywu/multiset/data/HCCDB/HCC/out_edgelist_Dataset1_HCCDB1.txt
/home/gywu/multiset/data/HCCDB/HCC/out_edgelist_Dataset2_HCCDB13.txt
/home/gywu/multiset/data/HCCDB/HCC/out_edgelist_Dataset3_HCCDB15.txt
/home/gywu/multiset/data/HCCDB/HCC/out_edgelist_Dataset4_HCCDB17.txt
/home/gywu/multiset/data/HCCDB/HCC/out_edgelist_Dataset5_HCCDB18.txt
/home/gywu/multiset/data/HCCDB/HCC/out_edgelist_Dataset6_HCCDB3.txt
/home/gywu/multiset/data/HCCDB/HCC/out_edgelist_Dataset7_HCCDB4.txt
/home/gywu/multiset/data/HCCDB/HCC/out_edgelist_Dataset8_HCCDB6.txt
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
