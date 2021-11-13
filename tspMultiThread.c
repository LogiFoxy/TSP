#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>

// N is number of total nodes on the graph or the cities in the map
#define N 100
// total number of threads
#define NUM_THREADS 2

// value for representing infinity, equals to max integer value
#define INF 2147483647

// helper structure to hold path
typedef struct pair {
  int first;
  int second;
} Pair;

// State Space Tree nodes
typedef struct node {

	// stores the reduced matrix
	int reducedMatrix[N][N];

	// stores the lower bound
	int cost;

	//stores current city number
	int vertex;

	// stores number of cities visited so far
	int level;

  int lastPosInPath;

  // stores the path from which we arrived to this node
  Pair path[N];
} Node;

// simulates a priority queue, where elements are accessed by their priority
// (in this case, lower priority comes first)
typedef struct pq_node {

	// Lower values indicate higher priority
	int priority;

  // pointer to the next node in the priority queue
	struct pq_node* next;

  // node data (node from graph) of this priority queue node
  struct node* node;

} PQ_node;

// helper structures to pass arguments to threads
typedef struct thread_args {
  // row of reduced matrix to be edited
  int i;

  // column of reduced matrix to be edited
  int j;

  // node with minimum cost extracted from the priority queue
  struct node* min;

} Thread_args;

// the priority queue is a shared structure between threads
// and hence it is accessed with a lock
PQ_node* pq;
pthread_mutex_t lock;

// definition of thread's run function
void *threadRun(void *args);

// function to create a new priority queue node
PQ_node* newPqNode(Node* node, int p) {
	PQ_node* temp = malloc(sizeof *temp);
	temp->node = node;
	temp->priority = p;
	temp->next = NULL;

	return temp;
}

// return the value at head of the priority queue
Node* peek(struct pq_node** head) {
	return (*head)->node;
}

// removes the element with the
// highest priority (=lowest value) form the list
void pop(struct pq_node** head) {
	struct pq_node* temp = *head;
	(*head) = (*head)->next;
	free(temp);
}

// function to push to the priority queue according to priority
void push(struct pq_node** head, Node* d, int p) {

  if((*head) == NULL ) {
    (*head) = newPqNode(d, p);
  }

	PQ_node* start = (*head);

	// create new Node
	PQ_node* temp = newPqNode(d, p);

  // insert new node before head
	if ((*head)->priority > p) {
		temp->next = *head;
		(*head) = temp;
	}
  // Traverse the list and find a position to insert new node
	else {
		while (start->next != NULL &&
			start->next->priority > p) {
			start = start->next;
		}
		// Either at the ends of the list
		// or at required position
		temp->next = start->next;
		start->next = temp;
	}
}

// function to check if priority queue is empty
int isEmpty(PQ_node** head) {
	return (*head) == NULL;
}

// function to create a new node, corresponds to visiting
// city in row from city in col
Node* newNode(int parentMatrix[N][N], Pair path[], int level, int row, int col) {

  Node* node = malloc(sizeof *node);
  node->lastPosInPath = 0;

	// stores ancestors edges of state space tree
  int pathSize = (int) (sizeof(path[0])/ sizeof(path[0]));
  for(int i=0;i<N;i++) {
    node->path[node->lastPosInPath].first = path[i].first;
    node->path[node->lastPosInPath].second = path[i].second;
    if(path[i].first!=-1 && path[i].second) node->lastPosInPath++;
  }

	// skip for root node
	if (level != 0) {
    // add current edge to path
    node->path[node->lastPosInPath].first=row;
    node->path[node->lastPosInPath].second=col;
    node->lastPosInPath++;
  }

	// copy data from parent node to current node
	memcpy(node->reducedMatrix, parentMatrix, sizeof node->reducedMatrix);

	// Change all entries of row row and column col to infinity
	// skip for root node
	for (int k = 0; level != 0 && k < N; k++) {
		// set outgoing edges for city row to infinity
		node->reducedMatrix[row][k] = INF;

		// set incoming edges to city col to infinity
		node->reducedMatrix[k][col] = INF;
	}

	// Set (col, 0) to infinity
	// here start node is 0
	node->reducedMatrix[col][0] = INF;

	// set number of cities visited so far
	node->level = level;

	// assign current city number
	node->vertex = col;

	return node;
}

// function to reduce each row by subtracting the minimum element
// from all elements in the row
void rowReduction(int reducedMatrix[N][N], int row[N]) {
	// initialize row array to INF
  for(int i = 0; i < N; i++) {
    row[i] = INF;
  }

	// row[i] contains minimum in row i of reducedMatrix
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			if (reducedMatrix[i][j] < row[i])
				row[i] = reducedMatrix[i][j];

	// reduce the minimum value from each element in each row
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			if (reducedMatrix[i][j] != INF && row[i] != INF)
				reducedMatrix[i][j] -= row[i];
}

// function to reduce each column by subtracting the minimum element
// from all elements in the column
void columnReduction(int reducedMatrix[N][N], int col[N]) {
	// initialize col array to INF
  for(int i = 0; i < N; i++) {
    col[i] = INF;
  }

	// col[j] contains minimum in col j
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			if (reducedMatrix[i][j] < col[j])
				col[j] = reducedMatrix[i][j];

	// reduce the minimum value from each element in each column
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			if (reducedMatrix[i][j] != INF && col[j] != INF)
				reducedMatrix[i][j] -= col[j];
}

// function to get the lower bound on
// on the path starting at current min node
int calculateCost(int reducedMatrix[N][N]) {
	// initialize cost to 0
	int cost = 0;

	// perform row reduction
	int row[N];
	rowReduction(reducedMatrix, row);

	// perform column
	int col[N];
	columnReduction(reducedMatrix, col);

	// the total expected cost
	// is the sum of all reductions
	for (int i = 0; i < N; i++)
		cost += (row[i] != INF) ? row[i] : 0,
			cost += (col[i] != INF) ? col[i] : 0;

	return cost;
}

// solve Traveling Salesman Problem using Branch and Bound
// using threads everytime it is needed to calculate the cost
// of a series of new nodes (expect for root node)
int solve(int costMatrix[N][N]) {

  struct pair path[N];
  for(int i=0; i<N; i++) {
    path[i].first = -1;
    path[i].second = -1;
  }

	// create a root node and calculate its cost
  Node* root = newNode(costMatrix, path, 0, -1, 0);

	// get the lower bound of the path starting at node 0
	root->cost = calculateCost(root->reducedMatrix);

  // create a priority queue to store live nodes of search tree
  pq = newPqNode(root, root->cost);

  // an array for the threads to be created
  pthread_t threads[NUM_THREADS];

	// find the live node with least cost, add its children to the priority queue
	// and finally delete it from the queue
	while (!isEmpty(&pq)) {
		// find the live node with least estimated cost
    Node *min = peek(&pq);

		// the found node is deleted from the queue of live nodes
    pop(&pq);
    // i stores current city number
		int i = min->vertex;

    //if all cities are visited
		if (min->level == N - 1) {
			// return to starting city
      min->path[min->lastPosInPath].first = i;
      min->path[min->lastPosInPath].second = 0;
      min->lastPosInPath++;

      // print list of cities visited;
      printf("Path to minimum cost:\n");
      for (int i = 0; i < min->lastPosInPath; i++) {
        printf("%d -> %d\n", min->path[i].first + 1, min->path[i].second + 1);
      }

			// return optimal cost
			return min->cost;
		}

    // do for each child of min
		// (i, j) forms an edge in space tree

    for (int j = 0; j < NUM_THREADS; j++) {
      Thread_args *args = malloc(sizeof *args);
      args->i = i;
      args->j = j;
      args->min = min;

      // create NUM_THREADS threads and run the threadRun function
      pthread_create(&threads[j], NULL, threadRun, args);

      // wait for threads to finish before proceeding
      pthread_join(threads[j], NULL);
    }

		// free node as we have already stored edges (i, j) in vector
		free(min);
	}
  return 0;
}

// the function to run when threads are created
// the function takes as parameters the node with minimum cost, extracted from
// the priority queue, i which represents the row to be edited from the costMatrix
// and its position in the thread array (j) so as to calculate the columns to edit
// and pushes child nodes (of the minimum cost node) to the priority queue
void *threadRun(void *args) {

  // extract arguments from  *args
  Thread_args *ta = (Thread_args *) args;
  int i = ta->i;
  int id = ta->j;
  // starting column for this thread
  int start = id*(N/NUM_THREADS);
  // final column for this thread
  int stop = start + (N/NUM_THREADS);
  if(id==NUM_THREADS-1) stop = N;

  for (int j = start; j < stop; j++) {
    if (ta->min->reducedMatrix[i][j] != INF) {
      // create a child node and calculate its cost
      Node* child = newNode(ta->min->reducedMatrix, ta->min->path,
        ta->min->level + 1, i, j);

      /* Cost of the child =
        cost of parent node +
        cost of the edge(i, j) +
        lower bound of the path starting at node j
      */
      child->cost = ta->min->cost + ta->min->reducedMatrix[i][j]
            + calculateCost(child->reducedMatrix);

      // mutual exclusion for shared structure
      pthread_mutex_lock(&lock);
      // Add child to list of live nodes
      push(&pq,child,child->cost);
      pthread_mutex_unlock(&lock);
    }
  }

  free(ta);
  return NULL;
}

// main function
int main() {

  // cost matrix for traveling salesman problem, representing the graph
  int costMatrix[N][N];

  // seed random function
  srand(time(NULL));

  // fill costMatrix with random values (values up to 100)
  for(int i=0; i<N; i++) {
    for(int j=0; j<N; j++) {
      if(i==j) {
          costMatrix[i][j] = INF;
      }
      else {
        costMatrix[i][j] = rand()%100;
      }
    }
  }

  // Values to test validity of program
  // !! if using this matrix, don't forget to change value of N to 5
  // int costMatrix[N][N] =
  // {
  // 	{ INF, 20, 30, 10, 11 },
  // 	{ 15, INF, 16, 4,	2 },
  // 	{ 3,	5, INF, 2,	4 },
  // 	{ 19,	6,	18,	INF, 3 },
  // 	{ 16,	4,	7,	16,	INF }
  // };

  clock_t begin = clock();
  int cost = solve(costMatrix);
  printf("Total cost is: %d\n",cost);
  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Total time: %f seconds\n", time_spent);

	return 0;
}
