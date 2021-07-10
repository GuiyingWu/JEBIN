#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include <iostream>


#define MAX_STRING 200
#define SIGMOID_BOUND 6

const int neg_table_size = 1e8;
const int hash_table_size = 30000000;
const int SIGMOID_TABLE_SIZE = 1000;

typedef float real;

typedef Eigen::Matrix< real, Eigen::Dynamic,
	Eigen::Dynamic, Eigen::RowMajor | Eigen::AutoAlign >
	BLPMatrix;

typedef Eigen::Matrix< real, 1, Eigen::Dynamic,
	Eigen::RowMajor | Eigen::AutoAlign >
	BLPVector;

struct struct_node {
	char *word;
};

struct hin_nb {
	int nb_id;
	double eg_wei;
	char eg_tp;
};

class sampler
{
	long long n;
	long long *alias;
	double *prob;

public:
	sampler();
	~sampler();

	void init(long long ndata, double *p);
	long long draw(double ran1, double ran2);
};

class line_node
{
protected:
	int node_max_size, vector_size;
	char node_file[MAX_STRING];
	int *node_hash;
	real *_vec;
	int get_hash(char *word);
	int add_node(char *word);
	
public:
	struct struct_node *node;
	int node_size;
	Eigen::Map<BLPMatrix> vec;

	line_node();
	~line_node();

	friend class line_hin;
	friend class line_trainer;

	void init_from_file(char *file_name, int vector_dim);
	void init_from_list(const std::vector<std::string> &node_list, int vector_dim);
	void init_vec_from_consensus(line_node *p_consensus_nodes);

	int search(char *word);
	void output(char *file_name, int binary);
};

class line_hin
{
protected:
	char hin_file[MAX_STRING];
	int vector_size;
	std::vector<hin_nb> *hin;
	std::vector<std::string> node_u_list;
	std::vector<std::string> node_v_list;
	long long hin_size;
	std::vector<std::string>::iterator it;

public:
	line_node *node_u, *node_v;
	
	line_hin();
	~line_hin();

	friend class line_trainer;

	void init(char *file_name, int vector_dim);
};

class line_trainer
{
protected:
	line_hin *phin;

	int *u_nb_cnt; int **u_nb_id; double **u_nb_wei;
	double *u_wei, *v_wei;
	sampler smp_u, *smp_u_nb;
	int neg_samples, *neg_table;

public:
	line_trainer();
	~line_trainer();

	void init(line_hin *p_hin, int negative);
	void train_sample(line_node *p_consensus_nodes, real rho, real gamma, real *_error_vec, real(*FastSigmoid)(real x), double(*func_rand_num)(), unsigned long long &rand_index);
};
