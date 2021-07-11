#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <math.h>
#include <time.h>
#include <sys/time.h> //gettimeofday
#include <pthread.h>
#include <gsl/gsl_rng.h>
#include <Eigen/Dense>
#include <iostream>
#include <vector>


#include "linelib.h"

char consensus_nodes_file[MAX_STRING], network_filenames_file[MAX_STRING], output_directory[MAX_STRING];
int num_network, binary = 0, num_threads = 1, vector_size = 200, negative = 5;
long long samples = 1, edge_count_actual;
real init_rho = 0.025, rho, gamma_value= 1.0, *sigmoid_table;


const gsl_rng_type * gsl_T;
gsl_rng * gsl_r;


std::vector<std::string> network_filenames_list;
line_node consensus_nodes;

std::vector<line_node> consensus_nodes_net_list;
std::vector<line_hin> network_in_list;
std::vector<line_trainer> trainer_list;


double func_rand_num()
{
	return gsl_rng_uniform(gsl_r);
}

/* Fastly compute sigmoid function */
void InitSigmoidTable()
{
	sigmoid_table = (real *)malloc((SIGMOID_TABLE_SIZE + 1) * sizeof(real));
	for (int i = 0; i < SIGMOID_TABLE_SIZE; i++)
	{
		sigmoid_table[i] = exp((i / (real)SIGMOID_TABLE_SIZE * 2 - 1) * SIGMOID_BOUND); // Precompute the exp() table
		sigmoid_table[i] = sigmoid_table[i] / (sigmoid_table[i] + 1);                   // Precompute f(x) = x / (x + 1)
	}
}

real FastSigmoid(real x)
{
	if (x > SIGMOID_BOUND) return 1;
	else if (x < -SIGMOID_BOUND) return 0;
	int i = (x + SIGMOID_BOUND) * SIGMOID_TABLE_SIZE / SIGMOID_BOUND / 2;
	return sigmoid_table[i];
}


void *TrainModelThread(void *id) 
{
	long long edge_count = 0, last_edge_count = 0;
	unsigned long long next_random = (long long)id;
	real *error_vec = (real *)calloc(vector_size, sizeof(real));
	int i;

	while (1)
	{
		if (edge_count > samples / num_threads + 2) break;

		if (edge_count - last_edge_count>10000)
		{
			edge_count_actual += edge_count - last_edge_count;
			last_edge_count = edge_count;
			printf("%crho: %f Progress: %.3lf%%", 13, rho, (real)edge_count_actual / (real)(samples + 1) * 100);
			fflush(stdout);
			rho = init_rho * (1 - edge_count_actual / (real)(samples + 1));
			if (rho < init_rho * 0.0001) rho = init_rho * 0.0001;
		}
		
		//gradient descent
		for(i=0; i<num_network; i++)
		{
			trainer_list[i].train_sample(&consensus_nodes, rho, gamma_value, error_vec, FastSigmoid, func_rand_num, next_random);
		}

		edge_count += num_network;	
	}
	free(error_vec);
	pthread_exit(NULL);
}

void TrainModel() {
	long a;
	int i;
	pthread_t *pt = (pthread_t *)malloc(num_threads * sizeof(pthread_t));
	rho = init_rho;
	
	network_in_list.resize(num_network);
	trainer_list.resize(num_network);

	consensus_nodes.init_from_file(consensus_nodes_file, vector_size);
	
	
	//initializing vector of u and v nodes for each network
	FILE *fi;
	char net_filenames[MAX_STRING];
	std::string str_net_filenames;
	fi = fopen(network_filenames_file, "rb");
	while (fscanf(fi, "%s\n", net_filenames) == 1)
	{
		printf("Network file: %s\n", net_filenames);
		str_net_filenames = net_filenames;
		network_filenames_list.push_back(str_net_filenames);
	}
	fclose(fi);
	
	for(i=0; i<num_network; i++)
	{
		str_net_filenames = network_filenames_list[i];
		strcpy(net_filenames, str_net_filenames.c_str());
		network_in_list[i].init(net_filenames, vector_size);
		//reinitializing vector of u nodes using consensus u nodes' vector for each network
		network_in_list[i].node_u->init_vec_from_consensus(&consensus_nodes);
		printf("net %d u Node size: %d\n", i, network_in_list[i].node_u->node_size);
		printf("net %d v Node size: %d\n", i, network_in_list[i].node_v->node_size);
	}
	
	
	//network structure processing ang initializing neg_table
	for(i=0; i<num_network; i++)
	{
		trainer_list[i].init(&network_in_list[i], negative);
	}
	
	InitSigmoidTable();
	
	
	//CPU time
	clock_t start = clock();
	//real time
	timeval t_start, t_end;
	gettimeofday(&t_start, NULL);
	
	printf("Training process:\n");
	for (a = 0; a < num_threads; a++) pthread_create(&pt[a], NULL, TrainModelThread, (void *)a);
	for (a = 0; a < num_threads; a++) pthread_join(pt[a], NULL);
	printf("\n");

	//real time
	gettimeofday(&t_end, NULL);
	double delta_t = (t_end.tv_sec-t_start.tv_sec) + (t_end.tv_usec-t_start.tv_usec)/1000000.0;
	printf("Total real time: %lf second\n", delta_t);
	
	//CPU time
	clock_t finish = clock();
	printf("Total CPU time: %lf second\n", (double)(finish - start) / CLOCKS_PER_SEC);
	printf("Number of threads: %d\n", num_threads);
	
	
	char temp_file[MAX_STRING];
	std::string output_dir = output_directory;

	//output consensus nodes' vector
	strcpy(temp_file, (output_dir+"output_u_consensus.txt").c_str());
	consensus_nodes.output(temp_file, binary);
	
	
	std::string str_i;
	char char_i[10];
	for(i=0; i<num_network; i++)
	{
		//output u nodes' vector
		sprintf(char_i, "%d", i+1);
		str_i = char_i;
		strcpy(temp_file, (output_dir + "output_u_net" + str_i + ".txt").c_str());
		network_in_list[i].node_u->output(temp_file, binary);

		//output v nodes' vector
		strcpy(temp_file, (output_dir + "output_v_net" + str_i + ".txt").c_str());
		network_in_list[i].node_v->output(temp_file, binary);
	}

}

int ArgPos(char *str, int argc, char **argv) {
	int a;
	for (a = 1; a < argc; a++) if (!strcmp(str, argv[a])) {
		if (a == argc - 1) {
			printf("Argument missing for %s\n", str);
			exit(1);
		}
		return a;
	}
	return -1;
}

int main(int argc, char **argv) {
	int i;
 	if (argc == 1) {
		printf("JEBIN: Joint Embedding of multiple BIpartite Networks\n\n");
		printf("Options:\n");
		printf("Parameters for training:\n");
		//input files
		printf("\t-consensus_nodes_file <file>\n");
		printf("\t\tThe consensus gene set\n");
		printf("\t-network_filenames_file <file>\n");
		printf("\t\tThe file contains all the bipartite network files' names (using absolute paths) \n");
		printf("\t-num_network <int>\n");
		printf("\t\tThe number of networks \n");
		//output directory
		printf("\t-output_directory <file>\n");
		printf("\t\tThe output directory of all kinds of the nodes' vectors \n");
		//parameters
		printf("\t-binary <int>\n");
		printf("\t\tSave the resulting vectors in binary mode; default is 0 (off) \n");
		printf("\t-size <int>\n");
		printf("\t\tThe dimension of the embedding vectors; default is 200 \n");
		printf("\t-negative <int>\n");
		printf("\t\tThe number of negative examples used in negative sampling; default is 5 \n");
		printf("\t-samples <float>\n");
		printf("\t\tThe total number of training samples (*Million) \n");
		printf("\t-threads <int>\n");
		printf("\t\tThe total number of threads used; the default is 1 \n");
		printf("\t-gamma <float>\n");
		printf("\t\tThe regularizing coefficient; the default is 1.0 \n");
		printf("\t-rho <float>\n");
		printf("\t\tThe starting value of the learning rate; the default is 0.025 \n");
		printf("\nExamples:\n");
		printf("./JEBIN_linux -consensus_nodes_file consensus_nodes_list.txt -network_filenames_file network_filenames_list.txt -num_network 2 -output_directory output/ -binary 0 -size 200 -negative 5 -samples 200 -threads 10 -gamma 1 -rho 0.025 \n\n");
		return 0;
	} 

	output_directory[0] = 0;

	if ((i = ArgPos((char *)"-consensus_nodes_file", argc, argv)) > 0) strcpy(consensus_nodes_file, argv[i + 1]);
	if ((i = ArgPos((char *)"-network_filenames_file", argc, argv)) > 0) strcpy(network_filenames_file, argv[i + 1]);
	if ((i = ArgPos((char *)"-num_network", argc, argv)) > 0) num_network = atoi(argv[i + 1]);
	if ((i = ArgPos((char *)"-output_directory", argc, argv)) > 0) strcpy(output_directory, argv[i + 1]);
	if ((i = ArgPos((char *)"-binary", argc, argv)) > 0) binary = atoi(argv[i + 1]);
	if ((i = ArgPos((char *)"-size", argc, argv)) > 0) vector_size = atoi(argv[i + 1]);
	if ((i = ArgPos((char *)"-negative", argc, argv)) > 0) negative = atoi(argv[i + 1]);
	if ((i = ArgPos((char *)"-samples", argc, argv)) > 0) samples = atof(argv[i + 1])*(long long)(1000000);
	if ((i = ArgPos((char *)"-gamma", argc, argv)) > 0) gamma_value = atof(argv[i + 1]);
	if ((i = ArgPos((char *)"-rho", argc, argv)) > 0) init_rho = atof(argv[i + 1]);
	if ((i = ArgPos((char *)"-threads", argc, argv)) > 0) num_threads = atoi(argv[i + 1]);
	
	printf("parameter: %s \n%s \n%s \n%d \n\n", consensus_nodes_file, network_filenames_file, output_directory, vector_size);
	
	gsl_rng_env_setup();
	gsl_T = gsl_rng_rand48;
	gsl_r = gsl_rng_alloc(gsl_T);
	gsl_rng_set(gsl_r, 314159265);

	TrainModel();
	return 0;
}