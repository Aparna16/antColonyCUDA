#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <iostream>
#include <curand.h>
#include <curand_kernel.h>

#define Infinity 65536 /* pow (2, 16) */

/* 
	randdouble()

	Retorna um numero (double) aleatorio entre 0.0f e 1.0f. 
	
	Parametros: 	
	Saida:
		numero aleaorio entre 0.0 e 1.0.
*/
#define randdouble() ((double)rand()/(double)RAND_MAX) 

/*
	randomize()

	Atualiza o gerador de numeros pseudo-aletatorios. 

	Parametros: 
	Saida:		
*/

#define randomize() srand((unsigned)time(NULL))

/* 
	index()

	Mapeia uma posicao de uma matriz (2D) para um indice de um vetor (1D). 
	
	Parametros:
		length: numero de colunas da matriz
		line: indice da linha 
		column: indice da coluna 

	Saída:
		indice mapeado
*/ 
#define index(length,line,column) (column + line * length) 

using namespace std;

const int NUMBER_OF_ITERATIONS = 50;
const double INIT_PHEROMONE_AMOUNT = 1.0;
const double EVAPORATION_RATE = 0.5;
const double ALFA = 1; /* Influencia da trilha de feromonios */
const double BETA = 2; /* Influencia da informacao heuristica */

/* 
	load_instance()

	Inicializa uma instancia ( numero de cidades e matriz de distancias ) do TSP. 

	Parametros:
		filename: nome do arquivo 
		n_cities: numero de cidades ( passagem por referencia ) 

	Saida: 
		matriz de distancias ( distancias euclidianas )
*/
int *load_instance ( char const *filename, int &n_cities ); 

/* 
	calculate_pathcost()

	Calcula o custo (soma dos custos de todas as arestas) de um caminho . 

	Parametros: 
		distances: matriz de distancias 
		path: caminho ( solucao ) 
		n_cities: numero de cidades 

	Saida:
		custo ( soma de todas as distancias ) do caminho
*/
int calculate_pathcost ( int *distances, int *path, int n_cities ); 

/* 
	best_solution()

	Retorna a melhor entre as solucoes geradas. 

	Parametros: 
		ants: matriz de solucoes 
		distances: matriz de distancias ( entre cidades ) 
		n_ants: numero de formigas 
		n_cities: numero de cidades

	Saida:
		melhor solucao encontrada 
*/
int *best_solution ( int *tours, int *distances, int n_ants, int n_cities ); 

/* 
	evaporate()

	Atualiza a matriz de feromonios aplicando evaporacao. 
	Para cada vertice, multiplica-se a taxa de evaporacao ( EVAPORATION_RATE ). 

	Parametros: 
		pheromones: matriz de feromonios 
		n_cities: numero de cidades 

	Saida:
		matriz de feromonios atualizada
*/
void evaporate ( double *pheromones, int n_cities ); 

/* 
	reinforce()

	Atualiza a matriz de feromonios. 
	Para cada vertice da melhor solucao corrente, adiciona-se uma quantidade de feromonios. 

	Parametros: 
		pheromones: matriz de feromonios 
		distances: matriz de distancias 
		min_path: caminho minimo ( melhor solucao ) encontrado 
		n_cities: numero de cidades 

	Saida: 
		matriz de feromonios atualizada
*/
void reinforce ( double *pheromones, int *distances, int *min_path, int n_cities ); 

/* 
	run()

	Executa o algoritmo de colonia de formigas 

	Parametros: 	
		distances: matriz de distancias 
		n_cities: numero de cidades 
		n_ants: numero de formigas 

	Saida: 
		melhor entre as solucoes encontradas por todas as formigas 
*/
int *run ( int *distances, int n_cities, int n_ants ); 

__global__ void cuda_evaporate ( double *pheromones, int n_cities, double evap_rate );

__global__ void cuda_reinforce ( double *pheromones, int *distances, int *path, int n_cities, double amount );

__global__ void cuda_construct_tour (int *tours, int *visited, double *choiceinfo, double *probs, int n_cities );

int main ( int argc, char *argv[] ) {
	randomize();

	char const *inputname, *outputname;

	if ( argc < 2 ) {
		cout << "Missing input arguments!" << endl;
		cout << "Program " << argv[0] << " takes exactly 3 arguments." << endl;
		return 1;
	}

	if ( argc > 3 ) {
		cout << "Too many arguments in program " << argv[0] << "!" << endl;
		cout << "It takes exactly 3 arguments." << endl;
		return 1;
	}

	cout << "Running " << argv[0] << " with arguments: ";

	for (int i = 1; i < argc; i++) 
		cout << argv[i] << " ";
	cout << endl;

	inputname = argv[1];

	if ( !argv[2] ) {
		outputname = "results/output.txt";
	} else {
		outputname = argv[2];
	}

	int n_cities; /* Numero de cidades */ 	
	int *distances;	/* Matriz de distancias (distancia euclidiana) */ 
	
	/*
		Inicializa a instancia. 
		Executa o algoritmo e calcula do custo da solucao. 
	*/
	distances 	= load_instance ( inputname, n_cities );
	int *solution = run ( distances, n_cities, 128 ); 
	int cost = calculate_pathcost ( distances, solution, n_cities ); 

	cout << "Writing results in file " << outputname << "!\n";
	ofstream output;
	output.open(outputname); 

	output << "Custo: " << cost << endl; 
	output << "Melhor solucao encontrada:\n"; 
	for(int i=0; i<n_cities; i++) 
		output << solution[i] << endl;
	
	cout << argv[0] << " exited with no errors."; 

	return 0; 
}

__global__ void cuda_evaporate ( double *pheromones, int n_cities, double evap_rate ) { 

	int edge_id = threadIdx.x + blockIdx.x * blockDim.x; 
	pheromones[ edge_id ] *= evap_rate;
} 

__global__ void cuda_reinforce ( double *pheromones, int *distances, int *path, int n_cities, double amount ) {

	int col_id = threadIdx.x + blockIdx.x * blockDim.x;

	int origin = path[col_id];
	int dest = path[col_id+1];

	pheromones[ index( n_cities, origin, dest ) ] += amount;
	pheromones[ index( n_cities, dest, origin ) ] += amount; 
}

__global__ void cuda_construct_tour (int *tours, int *visited, double *choiceinfo, double *probs, int n_cities ) {

	//int line_id = blockIdx.x; 
	int line_id = blockDim.x * blockIdx.x + threadIdx.x; 
	
	//extern __shared__ int shared_visited[]; 

	//for(int i = 0; i < n_cities; i++) {
		//shared_visited[ index ( n_cities, line_id, i ) ] = visited[ index ( n_cities, line_id, i ) ];
	//} 

	// __syncthreads();
	
	for (int step = 1; step < n_cities; step++) { 

		int current = tours[ index ( n_cities, line_id, step - 1 ) ];
		double sum_probs = 0.0;

		for(int i = 0; i < n_cities; i++) {
			if ( visited[ index ( n_cities, line_id, i) ] == 1 ) 
				probs[ index ( n_cities, line_id, i ) ] = 0.0;
			else {
				double current_prob = choiceinfo[ index( n_cities, current, i ) ]; 
				probs[ index ( n_cities, line_id, i ) ] = current_prob; 
				sum_probs += current_prob; 
			}
		}

		double random;
		curandState_t state;
		curand_init ( (unsigned long long) clock(), 0, 0, &state );
		random = curand_uniform ( &state ); 
		random *= sum_probs; 

		int next;
		double sum = probs[ index ( n_cities, line_id, 0 ) ];

		for(next = 0; sum < random; next++) {
			sum += probs[ index ( n_cities, line_id, next + 1 ) ];
		} 

		tours[ index ( n_cities, line_id, step ) ] = next; 
		visited[ index ( n_cities, line_id, next) ] = 1; 
	} 	
}

int *load_instance ( char const *filename, int &n_cities ) {
	cout << "Opening file " << filename << endl;

	ifstream tsp; 
	tsp.open (filename); 

	/*if ( ifstream == NULL ) {
		cout << "File " << filename << " not found!\n";
		exit(1);
	}*/
	
	tsp >> n_cities; 
	
	int* distances = (int *) malloc ( n_cities * n_cities * sizeof(int) ); 
	
	for (int i = 0; i < n_cities; i++) 
		for (int j = 0; j < n_cities; j++) 
			tsp >> distances[ index(n_cities, i, j) ]; 
	
	return distances; 
}

int calculate_pathcost ( int *distances, int *path, int n_cities ) {
	int cost = 0; 
	
	for (int i = 0; i < (n_cities - 1); i++) 
		cost += distances[ index(n_cities, path[i], path[i+1]) ]; 
	return cost; 
}

int *best_solution ( int *tours, int *distances, int n_ants, int n_cities ) {
	int *best_tour = &tours[0];
	for (int tour = 0; tour < n_ants; tour++) 
		if (calculate_pathcost(distances, &tours[index(n_cities, tour, 0)], n_cities) < calculate_pathcost(distances, best_tour, n_cities)) 
			best_tour = &tours[index(n_cities, tour, 0)];
	return best_tour;
}

void evaporate ( double *pheromones, int n_cities ) { 

	int size = n_cities * n_cities * sizeof(double); 
	double *pheromones_device; 
	cudaMalloc ( (void**) &pheromones_device, size); 
	
	cudaMemcpy (pheromones_device, pheromones, size, cudaMemcpyHostToDevice); 
	
	cuda_evaporate <<< n_cities, n_cities >>> ( pheromones_device, n_cities, EVAPORATION_RATE ); 
	
	cudaMemcpy (pheromones, pheromones_device, size, cudaMemcpyDeviceToHost); 
	
	cudaFree (pheromones_device); 
}

void reinforce ( double *pheromones, int *distances, int *path, int n_cities ) { 
	double amount = (double) ( 1.0f / (double) calculate_pathcost ( distances, path, n_cities ) );

	int size_path = n_cities * sizeof(int);
	int size_int = n_cities * n_cities * sizeof(int);
	int size_double = n_cities * n_cities * sizeof(double);

	int *path_device;
	int *distances_device;
	double *pheromones_device;
	
	cudaMalloc((void**)&path_device, size_path);
	cudaMalloc((void**)&distances_device, size_int);
	cudaMalloc((void**)&pheromones_device, size_double); 

	cudaMemcpy (path_device, path, size_path, cudaMemcpyHostToDevice);
	cudaMemcpy (distances_device, distances, size_int, cudaMemcpyHostToDevice); 
	cudaMemcpy (pheromones_device, pheromones, size_double, cudaMemcpyHostToDevice); 	

	cuda_reinforce <<< 1, n_cities - 1 >>> (pheromones_device, distances_device, path_device, n_cities, amount); 

	cudaMemcpy (distances, distances_device, size_int, cudaMemcpyDeviceToHost); 
	cudaMemcpy (pheromones, pheromones_device, size_double, cudaMemcpyDeviceToHost); 	
	
	cudaFree (path_device); 
	cudaFree (distances_device); 
	cudaFree (pheromones_device); 
}

int *run ( int *distances, int n_cities, int n_ants) { 

	int ph_size = n_cities * n_cities * sizeof(double);
	int tours_size = n_ants * n_cities * sizeof(int);
	int dist_size = n_cities * n_cities * sizeof(int);

	double *pheromones = (double*) malloc ( ph_size ); 
	int *tours = (int*) malloc ( tours_size ); /* Solucoes */
	int *visited = (int*) malloc ( tours_size ); /* Lista de cidades visitadas */
	double *choiceinfo = (double*) malloc ( ph_size ); 

	int *distances_device; /* Copia da GPU da matriz de distancias */
	int *tours_device; /* Copia da GPU da matriz de solucoes */
	int *visited_device; /* Copia da GPU da matriz de cidades visitadas */ 
	double *choiceinfo_device; /* Copia da GPU da matriz de probabilidades (numeraodor) */
	double *probs; /* Matriz de probabilidades */

	cudaMalloc ( (void**) &distances_device, dist_size ); 
	cudaMalloc ( (void**) &tours_device, tours_size ); 
	cudaMalloc ( (void**) &visited_device, tours_size ); 	
	cudaMalloc ( (void**) &choiceinfo_device, ph_size ); 
	cudaMalloc ( (void**) &probs, ph_size); 
		
	cudaMemcpy ( distances_device, distances, dist_size, cudaMemcpyHostToDevice );

	/*
		Instancia-se a matriz de feromonios. 
		Inicialmente, todas as arestas possuem a mesma quantidade de feromonios ( INIT_PHEROMONE_AMOUNT ). 
	*/ 
	for (int i = 0; i < n_cities; i++) 
		for (int j = 0; j < n_cities; j++) 
			pheromones[ index(n_cities, i, j) ] = INIT_PHEROMONE_AMOUNT; 

	for (int iteration = 0; iteration < NUMBER_OF_ITERATIONS; iteration++) { 
		/*
			Reseta todos os caminhos ao inicio de cada iteracao. 
			Inicialmente, todas as posicoes encontram-se no infinito. 
		*/

		for (int i = 0; i < n_ants; i++) 
			for (int j = 0; j < n_cities; j++) 
				tours[ index(n_cities, i, j) ] = Infinity; 

		for (int i = 0; i < n_ants; i++) 
			for (int j = 0; j < n_cities; j++)
				visited[ index(n_cities, i, j) ] = 0; 

		/*
			Calcula o numerador da funcao de probabilidade. 
			Em cada iteracao, este valor eh o mesmo para cada formiga, o que encoraja sua execucao aqui, 
			aumentando o desempenho do algoritmo. 
		*/
		
		for (int i = 0; i < n_cities; i++) {
			for (int j = 0; j < n_cities; j++) {
				double edge_pherom 	= pheromones[ index(n_cities, i, j) ]; 
				double edge_weight 	= distances[index(n_cities, i, j) ]; 
				double prob 		= 0.0f;
				if ( edge_weight != 0.0f ) {
					prob = pow ( edge_pherom, ALFA ) * pow ( (1/edge_weight), BETA ); 
				} else {
					prob = pow ( edge_pherom, ALFA ) * pow ( Infinity, BETA ); 
				} 
				choiceinfo[index(n_cities, i, j)] = prob; 
			} 
		} 		
		
		cudaMemcpy ( choiceinfo_device, choiceinfo, ph_size, cudaMemcpyHostToDevice ); 

		for (int ant = 0; ant < n_ants; ant++) { 
			int step = 0;

			/*
				Uma cidade inicial eh selecionada aleatoriamente. 
			*/
			int init = rand() % n_cities;

			/*
				Atualiza o tour ( para cada formiga ).
			*/
			tours [ index ( n_cities, ant, step ) ] = init;
			
			/*
				Atualiza a memoria da formiga. 	
			*/
			visited [ index ( n_cities, ant, init ) ] = 1; 
		} 

		cudaMemcpy ( visited_device, visited, tours_size, cudaMemcpyHostToDevice ); 
		cudaMemcpy ( tours_device, tours, tours_size, cudaMemcpyHostToDevice ); 

		int blockDim = 8;
		int antsPerBlock = n_ants / blockDim;
		// int sharedMemorySize = n_ants * n_cities * sizeof(int);
		
		cuda_construct_tour <<< blockDim, antsPerBlock/*, sharedMemorySize */>>> ( tours_device, visited_device, choiceinfo_device, probs, n_cities ); 

		cudaMemcpy ( tours, tours_device, tours_size, cudaMemcpyDeviceToHost ); 
		cudaMemcpy ( visited, visited_device, tours_size, cudaMemcpyDeviceToHost );	

		evaporate ( pheromones, n_cities ); 
		int *best = best_solution ( tours, distances, n_ants, n_cities ); 
		reinforce ( pheromones, distances, best, n_cities ); 
	} 

	cudaFree ( distances_device ); 
	cudaFree ( tours_device ); 
	cudaFree ( visited_device ); 
	cudaFree ( choiceinfo_device );	
	cudaFree ( probs ); 
	
	int *best = best_solution ( tours, distances, n_ants, n_cities ); 
	return best; 
} 