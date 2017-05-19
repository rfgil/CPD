#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <mpi.h>

//                                          AVL TREE
// *********************************************************************************************************************
// *********************************************************************************************************************

#define SMALLER 1
#define GREATER 2
#define EQUAL 3

typedef struct avl_node {
  void * item;
  int height;
  struct avl_node * left;
  struct avl_node * right;
} AVLNode;

typedef struct tree {
  int n_elements;
  int (*compareItems) (void *, void *);
  int (*compareIdItem) (int, void *);
  AVLNode * root;
} AVLTree;


AVLTree * newAvlTree(int (*compareItems)(void *, void *), int (*compareIdItem)(int, void *)){
  AVLTree * new_tree;

  new_tree = (AVLTree *) malloc(sizeof(AVLTree));
  new_tree->n_elements = 0;
  new_tree->compareItems = compareItems;
  new_tree->compareIdItem = compareIdItem;
  new_tree->root = NULL;

  return new_tree;
}

static int max(int a, int b){
  return a > b ? a : b;
}

static int getHeight(AVLNode * avl_node){
  return avl_node == NULL ? 0 : avl_node->height;
}

static AVLNode * rightRotate(AVLNode * y){
  AVLNode * x = y->left;
  AVLNode * T2 = x->right;

  // Perform rotation
  x->right = y;
  y->left = T2;

  // Update heights
  y->height = max(getHeight(y->left), getHeight(y->right))+1;
  x->height = max(getHeight(x->left), getHeight(x->right))+1;

  // Return new root
  return x;
}

static AVLNode * leftRotate(AVLNode * x){
  AVLNode * y = x->right;
  AVLNode * T2 = y->left;

  // Perform rotation
  y->left = x;
  x->right = T2;

  //  Update heights
  x->height = max(getHeight(x->left), getHeight(x->right))+1;
  y->height = max(getHeight(y->left), getHeight(y->right))+1;

  // Return new root
  return y;
}

static AVLNode * insertAVLNode(AVLNode * root, AVLNode * new_avl_node, int (*compareItems)(void *, void *)){
  int balance;

  if(root == NULL){
    return new_avl_node;
  }

  // Inserção de acordo com uma binary search tree (Nodes inferiores à esquerda e superiores à direita)
  if (compareItems(new_avl_node->item, root->item) == SMALLER) {
    root->left = insertAVLNode(root->left, new_avl_node, compareItems);
  } else {
    root->right = insertAVLNode(root->right, new_avl_node, compareItems);
  }

  // Calcula informações para confirmar se a árvore está balanceada
  root->height = max(getHeight(root->left), getHeight(root->right)) + 1;
  balance = getHeight(root->left) - getHeight(root->right);

  // Faz as rotações necessárias para balancear a árvore
  if (balance > 1 && compareItems(new_avl_node->item, root->left->item) == SMALLER){
    //Left Left
    return rightRotate(root);
  }
  if (balance < -1 && compareItems(new_avl_node->item, root->right->item) == GREATER) {
    // Right Right
    return leftRotate(root);
  }
  if (balance > 1 && compareItems(new_avl_node->item, root->left->item) == GREATER) {
    // Left Right
    root->left =  leftRotate(root->left);
    return rightRotate(root);
  }
  if (balance < -1 && compareItems(new_avl_node->item, root->right->item) == SMALLER ) {
    // Right Left
    root->right = rightRotate(root->right);
    return leftRotate(root);
  }

  // Não é necessário alterações: Devolve o root original
  return root;
}

void insertAvlTree(AVLTree * tree, void * item){
  AVLNode * new_avl_node;

  new_avl_node = (AVLNode *) malloc(sizeof(AVLNode));
  new_avl_node->item = item;
  new_avl_node->height = 1;
  new_avl_node->right = NULL;
  new_avl_node->left = NULL;

  if (tree->n_elements == 0){
    tree->root = new_avl_node;
  } else {
    tree->root = insertAVLNode(tree->root, new_avl_node, tree->compareItems);
  }

  tree->n_elements += 1;
}

void * findAvlTree(AVLTree * tree, int identifier){
  AVLNode * aux;

  if (tree == NULL){
    return NULL;
  }

  aux = tree->root;
  while (aux != NULL){
    switch ( tree->compareIdItem(identifier, aux->item) ) {
      case GREATER:
        aux = aux->right;
        break;

      case SMALLER:
        aux = aux->left;
        break;

      case EQUAL:
        return aux->item;
        break;

      default:
        return NULL;
        break;
    }
  }

  return NULL;
}

static void freeAVLNodes(AVLNode * avl_node, void (*freeItem)(void *)){
  if (avl_node != NULL){
    freeAVLNodes(avl_node->left, freeItem);
    freeAVLNodes(avl_node->right, freeItem);
    freeItem(avl_node->item);
    free(avl_node);
  }
}

void freeAvlTree(AVLTree * tree, void (*freeItem)(void *)){
  if (tree != NULL){
    freeAVLNodes(tree->root, freeItem);
    free(tree);
  }
}



// Estruturas e funções relacionadas com a AVLTree
// Estrtura que representa um célula (será inserida na AVL tree)
typedef struct cell{
 int x, y, z;
 int neighbors[6];
} Cell;

void printCell(void * cell){
   printf("%d %d %d\n", ((Cell *)cell)->x, ((Cell *)cell)->y, ((Cell *)cell)->z);
}

int compareItems(void * a, void * b){
 if ( ((Cell*)a)->z < ((Cell*)b)->z){
   return SMALLER;
 } else {
   return GREATER;
 }
}

int compareIdItem(int a, void * b){
 if (a == ((Cell*)b)->z ){
   return EQUAL;
 } else if (a < ((Cell*)b)->z ){
   return SMALLER;
 } else {
   return GREATER;
 }
}


// *********************************************************************************************************************
//                                            LISTA
// *********************************************************************************************************************
// Tags das mensagens enviadas entre processadores
#define INIT_TAG 1
#define BORDERS_TAG 2
#define PRINT_TAG 3

// Numero máximo de células que um processador receberá
// (corresponde aproximadamente a 100KB tendo em conta que um célula são 3 ints)
#define MAX_BUFFER_SIZE 8500

typedef struct list_node{
  int x, y, z;
  struct list_node * next;
} ListNode;

typedef struct list_head{
    int counter;
    ListNode * first_node;
} List;

List * newList(){
  List * new_list;

  new_list = (List *) malloc(sizeof(List));
  new_list->counter = 0;
  new_list->first_node = NULL;

  return new_list;
}

void insertList(List * list, int x, int y, int z){
  ListNode * new_list_node;

  new_list_node = (ListNode *) malloc(sizeof(ListNode));
  new_list_node->x = x;
  new_list_node->y = y;
  new_list_node->z = z;
  new_list_node->next = list->first_node;

  list->counter ++;
  list->first_node = new_list_node;
}

void freeListNodes(ListNode * list_node){
  if (list_node != NULL){
    freeListNodes(list_node->next);
    free(list_node);
  }
}

void resetList(List * list){
  if (list != NULL){
    freeListNodes(list->first_node);
    list->counter = 0;
    list->first_node = NULL;
  }
}

void freeList(List * list){
  if (list != NULL){
    freeListNodes(list->first_node);
    free(list);
  }
}

int serializeAndFreeList(ListNode * node, int * buffer, int index){
  if (node != NULL){

    index = serializeAndFreeList(node->next, buffer, index);

    buffer[index++] = node->x;
    buffer[index++] = node->y;
    buffer[index++] = node->z;

    free(node);
  }

  return index;
}

int * serializeList(List * list){
  int * buffer;
  ListNode * aux = list->first_node;

  buffer = (int *) malloc(3*list->counter*sizeof(int));

  int i = 0;
  while(aux != NULL){
    buffer[i++] = aux->x;
    buffer[i++] = aux->y;
    buffer[i++] = aux->z;

    aux = aux->next;
  }

  return buffer;
}

void sendList(List * list, int dest, int tag){
  int * serialized;

  serialized = (int *) malloc(3*list->counter*sizeof(int));
  serializeAndFreeList(list->first_node, serialized, 0);
  list->first_node = NULL;

  //serialized = serializeList(list);
  MPI_Send(serialized, 3*list->counter, MPI_INT, dest, tag, MPI_COMM_WORLD);
  free(serialized);
  //resetList(list);
}

int * sendListNoBlock(List * list, int dest, int tag, MPI_Request * request){
  int * serialized;
  serialized = serializeList(list);
  MPI_Isend(serialized, 3*list->counter, MPI_INT, dest, tag, MPI_COMM_WORLD, request);
  //free(serialized);
  resetList(list);
  return serialized;
}



void printAvlRecursive(AVLNode * root, void (*printItem)(void *)){
  if (root != NULL){
    printAvlRecursive(root->left, printItem);
    printItem(root->item);
    printAvlRecursive(root->right, printItem);
  }
}

void printAvlTree(AVLTree * tree, void (*printItem)(void *)){
  if(tree != NULL){
    printAvlRecursive(tree->root, printItem);
  }
}

void AVLTreeToListRecursive(AVLNode * root, List * buffer_list){
  Cell * cell;
  if (root != NULL){
    AVLTreeToListRecursive(root->left, buffer_list);

    cell = (Cell *) root->item;
    insertList(buffer_list, cell->x, cell->y, cell->z);

    AVLTreeToListRecursive(root->right, buffer_list);
  }
}

// *********************************************************************************************************************
//                                            MAIN FILE
// *********************************************************************************************************************

// Tamanho da aresta do cubo
int SIZE;

// Tamanho das partições pelos processadores
int TOP_SIZE, BOT_SIZE; // Tamanho das partições
int TOP_COUNT, BOT_COUNT; //Quatidade de partições

int N_PROCESSORS, PROCESSOR_ID;

// Estados possiveis de uma célula
#define ALIVE 1
#define DEAD 0

// Pontos de vista dos vizinhos relativamente a uma célula
#define TOP 0
#define BOTTOM 1
#define X_RIGHT 2
#define X_LEFT 3
#define Y_RIGHT 4
#define Y_LEFT 5

// Normaliza as coordenadas, ou seja garante que o vizinho de uma fronteira é a fronteira oposta (sides wrap around)
#define checkVectorLimits(c) ((c) == SIZE ? 0 : ((c) == -1 ? SIZE-1 : (c)))

//Converte o ponto de vista de uma célula para com a sua vizinha, no ponto de vista da vizinha para a primeira
//POV original    {TOP,    BOTTOM, X_RIGHT, X_LEFT,  Y_RIGHT Y_LEFT}
const int POV[] = {BOTTOM, TOP,    X_LEFT,  X_RIGHT, Y_LEFT, Y_RIGHT};

// Em função do ponto de vista obtem as coordenadas da célula vizinha
#define POV_Z(i, z) ((i) == TOP     ? z+1 : ((i) == BOTTOM ? z-1 : (z)))
#define POV_X(i, x) ((i) == X_RIGHT ? x+1 : ((i) == X_LEFT ? x-1 : (x)))
#define POV_Y(i, y) ((i) == Y_RIGHT ? y+1 : ((i) == Y_LEFT ? y-1 : (y)))

// Obtem as coordenadas da célula vizinha e normaliza-as
#define GET_X(i, x) checkVectorLimits(POV_X(i, x))
#define GET_Y(i, y) checkVectorLimits(POV_Y(i, y))
#define GET_Z(i, z) checkVectorLimits(POV_Z(i, z))

// Permite calcular o endereço do vetor que repsenta a matriz xy
#define GET_INDEX(x, y) (SIZE)*(x) + (y)
#define INDEX_TO_X(index) index / SIZE;
#define INDEX_TO_Y(index) index % SIZE;

//#define POV_INDEX(i, index) GET_INDEX(POV_X(i, index / SIZE), POV_X(i, index % SIZE))

// Devolve o processador que deve processar determinado index
#define GET_PROCESSOR(index) ((index) < (BOT_COUNT * BOT_SIZE) ? (index) / BOT_SIZE : (((index) - BOT_COUNT * BOT_SIZE) / TOP_SIZE) + BOT_COUNT)

// Devolve o primeiro index que o processador 'processor' deve processar
#define FIRST_INDEX(processor) ((processor) < BOT_COUNT ? BOT_SIZE * PROCESSOR_ID : BOT_SIZE*BOT_COUNT + (PROCESSOR_ID - BOT_COUNT)*TOP_SIZE)

//#define NORMALIZE_X(x) ((x) == SIZE ? - 1 : ((x) == -1 ? SIZE : (x)))
#define NORMALIZE_X(x) (PROCESSOR_ID == 0 && x == SIZE - 1 ? -1 : (PROCESSOR_ID == N_PROCESSORS - 1 && x == 0 ? SIZE : (x)))

// Normaliza/Desnormaliza um index para determinado processador
#define NORMALIZE_INDEX(index, processor) (index) - FIRST_INDEX(processor) + SIZE
#define DENORMALIZE_INDEX(normalized_index, processor) (normalized_index) + FIRST_INDEX(processor) - SIZE

#define PROCESSOR_CHUNK_SIZE(processor) ((processor) < BOT_COUNT ? BOT_SIZE : TOP_SIZE)


void setPartitionsSize(){
  int n_squared = SIZE * SIZE;

  TOP_SIZE = (n_squared + (N_PROCESSORS - 1)) / N_PROCESSORS; // Equivalente a arrendondar para cima n_squared / n_processors
  BOT_SIZE = n_squared / N_PROCESSORS; // Arredondar para baixo é apenas truncar o resultado

  //if (!PROCESSOR_ID){
  //  printf("TOP_SIZE: %d\nBOT_SIZE: %d\n", TOP_SIZE, BOT_SIZE);
  //}

  if (TOP_SIZE == BOT_SIZE){
    TOP_COUNT = 0;
    BOT_COUNT = n_squared / BOT_SIZE;

    //if (!PROCESSOR_ID){
    //  printf("%d chunks de tamanho %d\n", BOT_COUNT, BOT_SIZE);
    //}
    return;
  }

  if(N_PROCESSORS > SIZE){
    // A cada processador com rank com id menor que size é atribuido um plano de z
    N_PROCESSORS = SIZE;
    BOT_COUNT = SIZE;
    BOT_SIZE = SIZE;

    //if (!PROCESSOR_ID){
    //  printf("%d chunks de tamanho %d\n", BOT_COUNT, BOT_SIZE);
    //}

    return;
  }

  for (TOP_COUNT = 1; TOP_COUNT < n_squared/BOT_SIZE; TOP_COUNT++){
    for (BOT_COUNT = 1; BOT_COUNT < n_squared/BOT_SIZE; BOT_COUNT++){
      if (TOP_SIZE*TOP_COUNT + BOT_SIZE*BOT_COUNT == n_squared){
        //if (!PROCESSOR_ID){
        //  printf("%d chunks de tamanho %d e %d chunks de tamanho %d\n", BOT_COUNT, BOT_SIZE, TOP_COUNT, TOP_SIZE);
        //}
        return;
      }
    }
  }

}

void insertNewCellInAVL(AVLTree * tree, int x, int y, int z){
  Cell * new_cell;

  new_cell = (Cell *) malloc(sizeof(Cell));

  new_cell->x = x;
  new_cell->y = y;
  new_cell->z = z;
  memset(new_cell->neighbors, -1, 6*sizeof(int)); //Inicializa a memória com -1

  insertAvlTree(tree, new_cell);
}

AVLTree ** loadMap(char * file_name){
  FILE * input_file;
  AVLTree ** map; // mapa para o processador 0
  List ** buffer_list;
  int x, y, z;

  int index, processor;

  input_file = fopen(file_name, "r");

  if (input_file == NULL || fscanf(input_file, "%d", &SIZE) != 1){
    SIZE = 0;
    MPI_Bcast(&SIZE, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Finalize();
    exit(0);
  }

  MPI_Bcast(&SIZE, 1, MPI_INT, 0, MPI_COMM_WORLD);
  setPartitionsSize();

  // Alocar estrutura para o processador de rank 0
  map = (AVLTree **) calloc(BOT_SIZE + 2*SIZE, sizeof(AVLTree *));
  buffer_list = (List **) malloc(N_PROCESSORS * sizeof(List *));

  for (index = 1; index<N_PROCESSORS; index++){ //Não é necessário inicializar para o proc 0
    buffer_list[index] = newList();
  }

  while(fscanf(input_file, "%d %d %d", &x, &y, &z) == 3){
    index = GET_INDEX(x, y);

    if (index < BOT_SIZE){
      // Esta célula pertence ao processador 0 (será imediamente adiconada ao mapa local)

      if (map[index + SIZE] == NULL){
        // No processador 0 normalizar o index corresponde a somar size
        map[index + SIZE] = newAvlTree(compareItems, compareIdItem);
      }

      insertNewCellInAVL(map[index + SIZE], x, y, z);

    } else {
      processor = GET_PROCESSOR(index);
      insertList(buffer_list[processor], x, y, z);

      if (buffer_list[processor]->counter == MAX_BUFFER_SIZE){
        sendList(buffer_list[processor], processor, INIT_TAG);
      }
    }
  }

  fclose(input_file);

  for (index = 1; index<N_PROCESSORS; index++){ //Não é necessário inicializar para o proc 0
    sendList(buffer_list[index], index, INIT_TAG);
  }

  return map;
}

AVLTree ** receiveMap(){
  AVLTree ** map;
  MPI_Status status;
  int * buffer;

  int index, counter, i;
  int x, y, z;

  MPI_Bcast(&SIZE, 1, MPI_INT, 0, MPI_COMM_WORLD);
  setPartitionsSize();

  if(SIZE == 0 || PROCESSOR_ID >= N_PROCESSORS){
    // Ocorreu um erro ao ler o ficheiro de entrada no processador 0
    MPI_Finalize();
    exit(0);
  }

  map = (AVLTree **) calloc(PROCESSOR_CHUNK_SIZE(PROCESSOR_ID) + 2*SIZE, sizeof(AVLTree *));

  buffer = (int*) malloc(3*MAX_BUFFER_SIZE*sizeof(int));

  do {
      MPI_Recv(buffer, MAX_BUFFER_SIZE, MPI_INT, 0, INIT_TAG, MPI_COMM_WORLD, &status);
      MPI_Get_count(&status, MPI_INT, &counter);

      i = 0;
      while(i < counter){
        x = buffer[i++];
        y = buffer[i++];
        z = buffer[i++];

        index = NORMALIZE_INDEX( GET_INDEX(x, y), PROCESSOR_ID );
        if (map[index] == NULL){ // Inicializa a AVL
          map[index] = newAvlTree(compareItems, compareIdItem);
        }

        insertNewCellInAVL(map[index], x, y, z);
      }

  } while(counter == MAX_BUFFER_SIZE);

  return map;
}

void receiveBorders(AVLTree ** map, int proc){
  MPI_Status status;
  int * buffer;
  int counter, i;
  int norm_x;
  int x, y, z;
  int norm_index;

  buffer = (int *) malloc(3*SIZE*SIZE*sizeof(int));

  MPI_Recv(buffer, 3*SIZE*SIZE, MPI_INT, proc, BORDERS_TAG, MPI_COMM_WORLD, &status);
  MPI_Get_count(&status, MPI_INT, &counter);

  i = 0;
  while(i < counter){
    x = buffer[i++];
    y = buffer[i++];
    z = buffer[i++];

    // Normaliza a coordenada x
    /*
    if (x == SIZE - 1){
      norm_x = -1;
    } else if(x == 0){
      norm_x = SIZE;
    } else {
      norm_x = x;
    }
    */

    norm_x = NORMALIZE_X(x);

    norm_index = NORMALIZE_INDEX(GET_INDEX(norm_x, y), PROCESSOR_ID);

    if (map[norm_index] == NULL){ // Inicializa a AVL
      map[norm_index] = newAvlTree(compareItems, compareIdItem);
    }

    // Não ha problema de x "não estar correcto", desde que a célula seja inserida na posição xy normalizada certa
    insertNewCellInAVL(map[norm_index], x, y, z);
  }
}

void sendBorders(AVLTree ** map){
  MPI_Request request;
  MPI_Status status;

  List ** buffer_list;
  int norm_index;
  int prev_processor, next_processor;

  int * send_buffer;

  // inicializa buffer
  buffer_list = (List **) malloc(2* sizeof(List *));

  buffer_list[0] = newList();
  buffer_list[1] = newList();

  prev_processor = PROCESSOR_ID == 0 ? N_PROCESSORS - 1 : PROCESSOR_ID - 1;
  next_processor = PROCESSOR_ID == N_PROCESSORS - 1 ? 0 : PROCESSOR_ID + 1;

  // Fronteira com processador anterior
  for(norm_index = SIZE; norm_index < SIZE + SIZE; norm_index ++){
    AVLTreeToListRecursive(map[norm_index] == NULL ? NULL : map[norm_index]->root, buffer_list[0]);
  }

  send_buffer = sendListNoBlock(buffer_list[0], prev_processor, BORDERS_TAG, &request);
  // Aproveita latência da transferencia

  // Fronteira com vizinho seguinte
  for(norm_index = PROCESSOR_CHUNK_SIZE(PROCESSOR_ID); norm_index < PROCESSOR_CHUNK_SIZE(PROCESSOR_ID) + SIZE; norm_index ++){
    AVLTreeToListRecursive(map[norm_index] == NULL ? NULL : map[norm_index]->root, buffer_list[1]);
  }

  receiveBorders(map, next_processor);
  MPI_Wait(&request, &status);
  free(send_buffer);

  send_buffer = sendListNoBlock(buffer_list[1], next_processor, BORDERS_TAG, &request);

  receiveBorders(map, prev_processor);
  MPI_Wait(&request, &status);
  free(send_buffer);

  freeList(buffer_list[0]);
  freeList(buffer_list[1]);
  free(buffer_list);
}

void addToNewGeneration(AVLTree ** next_map, int x, int y, int z, int nAliveNeighbors, int cell_state){
  if (nAliveNeighbors == 2 || nAliveNeighbors == 3 || (nAliveNeighbors == 4 && cell_state == ALIVE)){ // Condição para a proxima celula estar viva

    // Verifica se árvore corresponde às coordenadas x e y da celula já está criada
    if (next_map[NORMALIZE_INDEX(GET_INDEX(x,y), PROCESSOR_ID)] == NULL){
        next_map[NORMALIZE_INDEX(GET_INDEX(x,y), PROCESSOR_ID)] = newAvlTree(compareItems, compareIdItem); //Inicializa a AVL
    }

    // Aloca nova célula
    insertNewCellInAVL(next_map[NORMALIZE_INDEX(GET_INDEX(x,y), PROCESSOR_ID)], x, y, z);
  }
}

int visitNeighborsDead(AVLTree ** map, int x, int y, int z, int pov_caller){
  int i;
  int count = 1; //Uma célula morta só é chamada por uma célula viva!
  Cell *  neighbor_cell;

  int neighbor_x, neighbor_y, neighbor_z;
  int neighbor_index, index;

  index = NORMALIZE_INDEX(GET_INDEX(x,y), PROCESSOR_ID);

  if (index < SIZE || index >= SIZE + PROCESSOR_CHUNK_SIZE(PROCESSOR_ID)){
      // Se estiver fora dos indices de interesse do processador
      return 0; // É considerada como morta
  }

  for(i=0; i<6; i++){ // Visita os 6 vizinhos de uma célula
    if (i != pov_caller){
      neighbor_x = POV_X(i, x);
      neighbor_y = GET_Y(i, y);
      neighbor_z = GET_Z(i, z);

      neighbor_index = NORMALIZE_INDEX(GET_INDEX(neighbor_x,neighbor_y), PROCESSOR_ID);

      if (neighbor_index < 0 || neighbor_index > 2*SIZE + PROCESSOR_CHUNK_SIZE(PROCESSOR_ID) - 1){
          // Se estiver fora dos indices de interesse do processador
          continue; // É considerada como morta
      }

      neighbor_cell = (Cell *) findAvlTree(map[neighbor_index], neighbor_z);

      if (neighbor_cell != NULL){
        //printf("NEIGHBOR: %d %d %d\n", neighbor_cell->x, neighbor_cell->y, neighbor_cell->z);

        // A célula vizinha está viva!!
        // Informa o vizinho sobre o estado da célula que o chamou
        neighbor_cell->neighbors[POV[i]] = DEAD;
        count ++;
      }
    }
  }

  return count;
}

int visitNeighborsAlive(AVLTree ** map, AVLTree ** next_map, Cell * current_cell, int isCurrentCellFromBorder){
  int i, nAliveNeighbors;
  int count = 0;

  Cell *  neighbor_cell;
  int neighbor_x, neighbor_y, neighbor_z, neighbor_index;

  for(i=0; i<6; i++){ //Visita os 6 vizinhos de uma célula
    if (current_cell->neighbors[i] == -1){
      neighbor_x = POV_X(i, current_cell->x);
      neighbor_y = GET_Y(i, current_cell->y);
      neighbor_z = GET_Z(i, current_cell->z);

      neighbor_index = NORMALIZE_INDEX(GET_INDEX(neighbor_x, neighbor_y), PROCESSOR_ID);
      if (isCurrentCellFromBorder && (neighbor_index < SIZE || neighbor_index >= SIZE + PROCESSOR_CHUNK_SIZE(PROCESSOR_ID) )){
            // Se estiver fora dos indices de interesse do processador
            continue; // É considerada como morta
      }


      //printf("(%d) ALIVE: %d %d %d CALLER: %d %d %d\n", PROCESSOR_ID, neighbor_x, neighbor_y, neighbor_z, current_cell->x, current_cell->y, current_cell->z);

      neighbor_cell = (Cell *) findAvlTree(map[neighbor_index], neighbor_z);

      if (neighbor_cell != NULL){
        // A célula vizinha está viva!!
        // Informa o vizinho sobre o estado da célula que o chamou
        neighbor_cell->neighbors[POV[i]] = ALIVE;
        count ++;

      } else {
        //Avalia a célula morta (Todas as celulas mortas que renascem têm pelo menos uma célula viva como vizinha)
        current_cell->neighbors[i] = DEAD;

        nAliveNeighbors = visitNeighborsDead(map, neighbor_x, neighbor_y, neighbor_z, POV[i]);
        addToNewGeneration(next_map, neighbor_x, neighbor_y, neighbor_z, nAliveNeighbors, DEAD);
      }
    } else count += current_cell->neighbors[i];
  }

  return count;
}

void searchTree(AVLTree ** map, AVLTree ** next_map, AVLNode * root, int isCurrentCellFromBorder){
  Cell * current_cell;
  int nAliveNeighbors;

  if (root != NULL){
    // Percorre a AVL Tree
    searchTree(map, next_map, root->left, isCurrentCellFromBorder);
    searchTree(map, next_map, root->right, isCurrentCellFromBorder);

    current_cell = (Cell *) root->item;
    nAliveNeighbors = visitNeighborsAlive(map, next_map, current_cell, isCurrentCellFromBorder);
    addToNewGeneration(next_map, current_cell->x, current_cell->y, current_cell->z, nAliveNeighbors, ALIVE);
  }
}

void getNextGeneration(AVLTree ** map, AVLTree ** next_map){
  int xy;

  memset(next_map, 0, (PROCESSOR_CHUNK_SIZE(PROCESSOR_ID) + 2*SIZE)*sizeof(AVLTree *)); //Inicializa o vetor next_map com o valor NULL

  //for(xy= SIZE; xy<SIZE+PROCESSOR_CHUNK_SIZE(PROCESSOR_ID); xy++){
  for(xy= 0; xy<2*SIZE+PROCESSOR_CHUNK_SIZE(PROCESSOR_ID); xy++){
    if (map[xy] != NULL){
      searchTree(map, next_map, map[xy]->root, (xy < SIZE || xy >= SIZE + PROCESSOR_CHUNK_SIZE(PROCESSOR_ID) ));
    }
  }
}

void freeMapAVLTrees(AVLTree ** map){
  int xy;

  for(xy=0; xy<(PROCESSOR_CHUNK_SIZE(PROCESSOR_ID) + 2*SIZE); xy++){
      freeAvlTree(map[xy], free);
  }
}

void writeOutput(AVLTree ** map){
  MPI_Status status;
  int * buffer;
  int xy, i, x, y, z, counter;

  // Imprime o próprio conteudo em rank0
  for(xy=SIZE; xy<SIZE + PROCESSOR_CHUNK_SIZE(PROCESSOR_ID); xy++){
    printAvlTree(map[xy], printCell);
  }

  buffer = (int *) malloc(3*MAX_BUFFER_SIZE*sizeof(int));

  for (int proc = 1; proc<N_PROCESSORS; proc++){
    do {
      MPI_Recv(buffer, 3*MAX_BUFFER_SIZE, MPI_INT, proc, PRINT_TAG, MPI_COMM_WORLD, &status);
      MPI_Get_count(&status, MPI_INT, &counter);

      i = 0;
      while(i < counter){
        x = buffer[i++];
        y = buffer[i++];
        z = buffer[i++];

        printf("%d %d %d\n", x, y, z);
      }
    } while(counter == 3*MAX_BUFFER_SIZE);
  }

  free(buffer);
}

void sendOutput(AVLTree ** map){
  List * buffer_list;

  buffer_list = newList();


  for(int xy=SIZE; xy<SIZE + PROCESSOR_CHUNK_SIZE(PROCESSOR_ID); xy++){
    if (map[xy] != NULL){
      AVLTreeToListRecursive(map[xy]->root, buffer_list);

      if (buffer_list->counter == MAX_BUFFER_SIZE){
        sendList(buffer_list, 0, PRINT_TAG);
      }
    }
  }

  sendList(buffer_list, 0, PRINT_TAG);
}

int main(int argc, char *argv[]){
  int n_generations, i;
  AVLTree ** map, ** next_map,  ** aux_pointer;

  MPI_Init (&argc, &argv);

  MPI_Comm_rank (MPI_COMM_WORLD, &PROCESSOR_ID);
  MPI_Comm_size (MPI_COMM_WORLD, &N_PROCESSORS);


  if (argc < 3 || sscanf(argv[2],"%d", &n_generations) != 1){
    // printf("Número argumentos insuficiente\n");
    // printf("O segundo argumento não é um inteiro\n");
    MPI_Finalize();
    exit(0);
  }

  //n_generations = 1;

  MPI_Barrier(MPI_COMM_WORLD);

  if (!PROCESSOR_ID){
    // O processador 0 lê o mapa e distribui a informação
    map = loadMap(argv[1]);
    //map = loadMap("s5e50.in");
  } else {
    // Processadores restantes recebem informação de do processador 0
    map = receiveMap();
  }

  MPI_Barrier(MPI_COMM_WORLD);

  next_map = (AVLTree **) malloc((PROCESSOR_CHUNK_SIZE(PROCESSOR_ID) + 2*SIZE)*sizeof(AVLTree *));

  for (i=0; i<n_generations; i++){

    // Troca fronteiras
    sendBorders(map);

    // Calcula próxima geração
    getNextGeneration(map, next_map);

    // Troca ponteiros
    aux_pointer = map;
    map = next_map;
    next_map = aux_pointer;

    // Apaga fronteira anterior
    memset(map, 0, SIZE);
    memset(map + SIZE + PROCESSOR_CHUNK_SIZE(PROCESSOR_ID), 0, SIZE);

    /*
    if (PROCESSOR_ID){
      for(int i=0; i<(2*SIZE + PROCESSOR_CHUNK_SIZE(PROCESSOR_ID)); i++){
        if(map[i] != NULL){
          //printf("RANK%d index: %d\n", PROCESSOR_ID, i);

          printAvlTree(map[i], printCell);
          printf("%d\n", i);
        }
      }
      printf("--------------------------------\n");
    }
    */

    freeMapAVLTrees(next_map);
  }

  MPI_Barrier(MPI_COMM_WORLD);


  if (!PROCESSOR_ID){
    // print result
    writeOutput(map);
  } else {
    // send result
    sendOutput(map);
  }

  MPI_Finalize();

  //strcpy(file_name + strlen(file_name) - 2, "out\0"); //Substitui .in por .out

  //writeOutput(map);

  //freeMapAVLTrees(map);
  //free(map);
  //free(next_map);
  //free(file_name);

  return 0;
}
