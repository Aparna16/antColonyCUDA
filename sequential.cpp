#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <iostream>

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
	construct_tour()

	Constroi um caminho ( solucao ) viavel. 
	Para cada formiga, inicia-se um caminho a partir de uma cidade. 
	Para cada cidade da vizinhanca, eh calculada a probabilidade. 
	A proxima cidade a ser visitada eh escolhida probabilisticamente. 

	Parametros: 
		distances: matriz de distancias 
		pheromones: matriz de feromonios 
		choiceinfo: matriz de probibilidades 
		tour: tour ( solucao ) atual 
		visited: vetor binario de tamanho n_cities ( numero de cidades ) que indica o estado de casa cidade na iteracao
				 - 0: nao visitada 
				 - 1: visitada 
		step: indice da posicao atual do caminho ( tour ) percorrido ate o momento pela formiga 
		n_cities: numero de cidades 

	Saida: 
		caminho ( tour ) construido
*/
void construct_tour ( int *distances, double *pheromones, double *choiceinfo, int *tour, int *visited, int step, int n_cities ); 

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
int *run ( int *distances, int n_cities, int n_ants);

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

	int n_cities; 	/* Numero de cidades */ 	
	int *distances;	/* Matriz de distancias (distancia euclidiana) */ 
	
	/*
		Inicializa a instancia. 
		Executa o algoritmo e calcula do custo da solucao. 
	*/
	distances 	= load_instance ( inputname, n_cities );
	int *solution = run ( distances, n_cities, n_cities ); 
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

int *load_instance ( char const *filename, int &n_cities ) {
	cout << "Opening file " << filename << endl;

	ifstream tsp; 
	tsp.open (filename); 

	if(tsp == NULL) {
		cout << "File " << filename << " not found!\n";
		exit(1);
	}
	
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
	for (int i = 0; i < n_cities; i++) 
		for (int j = 0; j < n_cities; j++) 
			pheromones[ index(n_cities, i, j) ] *= EVAPORATION_RATE; 
}

void reinforce ( double *pheromones, int *distances, int *min_path, int n_cities ) { 
	double amount = (double) ( 1.0f / (double) calculate_pathcost ( distances, min_path, n_cities ) );
	for (int i = 0; i < (n_cities - 1); i++) { 
		int origin 	= min_path[i];
		int dest 	= min_path[i+1];
		pheromones[ index( n_cities, origin, dest ) ] += amount;
		pheromones[ index( n_cities, dest, origin ) ] += amount;
	}
}

void construct_tour ( int *distances, double *pheromones, double *choiceinfo, int *tour, int *visited, int step, int n_cities ) { 

	/* 
		Cidade atual. 
	*/
	int current = tour[step-1]; 
	
	/* 
		Vetor de probabilidades da vizinhanca. 
	*/
	double neighbour_probs[ n_cities ]; 
	
	/*
		Soma das probabilidades de todas as cidades da vizinhanca. 
	*/
	double sum_probs = 0.0f; 

	/*
		Para cada cidade viavel, eh verificada a probabilidade de esta ser visitada partindo-se da cidade atual. 
	*/
	for (int i = 0; i < n_cities; i++) {
		/*
			Se a cidade ja foi visitada, atribui-se probabilidade 0. 
		*/
		if ( visited[i] == 1 ) 
			neighbour_probs[i] = 0;
		else {
			double current_prob = choiceinfo[ index(n_cities, current, i) ]; 
			neighbour_probs[i] 	= current_prob; 
			sum_probs += current_prob; 
		}
	}

	/*
		Seleciona um numero aleatorio entre 0.0f e sum_probs. 
	*/
	double random = randdouble() * sum_probs;
	
	/*
		Indice da proxima a cidade a ser visitada
	*/
	int next; 
	double sum = neighbour_probs[0]; 
	/*
		Aplica a selecao por roulette wheel. 
	*/
	for (next = 0; sum < random; next++) {
		sum += neighbour_probs[next + 1];
	} 
	
	/*
		Adiciona a proxima cidade na solucao atual e atualiza a lista de cidades visitadas. 
	*/
	
	tour[step] = next; 
	visited[next] = 1; 
}

int *run ( int *distances, int n_cities, int n_ants) {
	/*
		Instancia-se a matriz de feromonios. 
		Inicialmente, todas as arestas possuem a mesma quantidade de feromonios ( INIT_PHEROMONE_AMOUNT ). 
	*/
	double pheromones[ n_cities * n_cities ];

	for (int i = 0; i < n_cities; i++) 
		for (int j = 0; j < n_cities; j++) 
			pheromones[ index(n_cities, i, j) ] = INIT_PHEROMONE_AMOUNT; 

	/*
		Instancia-se a matriz de solucoes. 
	*/
	int * tours = (int*) malloc( n_ants * n_cities * sizeof(int) );

	for (int iteration = 0; iteration < NUMBER_OF_ITERATIONS; iteration++) { 
		/*
			Reseta todos os caminhos ao inicio de cada iteracao. 
			Inicialmente, todas as posicoes encontram-se no infinito. 
		*/

		for (int i = 0; i <  n_ants; i++) 
			for (int j = 0; j < n_cities; j++) 
				tours[ index(n_cities, i, j) ] = Infinity; 

		/*
			Calcula o numerador da funcao de probabilidade. 
			Em cada iteracao, este valor eh o mesmo para cada formiga, o que encoraja sua execucao aqui, 
			aumentando o desempenho do algoritmo. 
		*/
		double choiceinfo [ n_cities * n_cities ]; 
		for (int i = 0; i < n_cities; i++) 
			for (int j = 0; j < n_cities; j++) {
				double edge_pherom 	= pheromones[ index(n_cities, i, j) ]; 
				double edge_weight 	= distances[index(n_cities, i, j) ]; 
				double prob 		= 0.0f;
				if ( edge_weight != 0.0f ) {
					prob = pow ( edge_pherom, ALFA ) * pow ( (1/edge_weight), BETA ); 
				} else {
					prob = pow ( edge_pherom, ALFA ) * pow ( Infinity, BETA ); 
					//prob = 1.0f;
				} 
				choiceinfo[index(n_cities, i, j)] = prob; 
			} 
		
		for (int ant = 0; ant < n_ants; ant++) { 
			int step = 0;
			/* 
				Informa as cidades ja visitadas pela formiga. 
				Para cada cidade, o valor eh 1 se a cidade ja foi visitada, e 0, caso contrario. 
				Incialmente, todas as posicoes sao 0. 
			*/
			int visited[ n_cities ]; 
			for (int i = 0; i < n_cities; i++) 
				visited[i] = 0;
			
			/*
				Caminho da formiga atual. 
				Sera construido ao longo da construcao do tour. 
			*/
			int *tour = &tours[ index( n_cities, ant, 0 ) ]; 

			/*
				Uma cidade inicial eh selecionada aleatoriamente. 
			*/
			int init = rand() % n_cities;

			tour[step] = init;
			visited[init] = 1;		

			for (step = 1; step < n_cities; step++) 
				construct_tour ( distances, pheromones, choiceinfo, tour, visited, step, n_cities ); 
		} 
		evaporate ( pheromones, n_cities );
		reinforce ( pheromones, distances, best_solution ( tours, distances, n_ants, n_cities ), n_cities );
	} 

	return best_solution ( tours, distances, n_ants, n_cities ); 
}