#include <stdlib.h>
#include <stdio.h>
#include <string.h>

//                                          AVL TREE
// *********************************************************************************************************************
// *********************************************************************************************************************

#define SMALLER 1
#define GREATER 2
#define EQUAL 3

typedef struct node {
  void * item;
  int height;
  struct node * left;
  struct node * right;
} Node;

typedef struct tree {
  int n_elements;
  int (*compareItems) (void *, void *);
  int (*compareIdItem) (int, void *);
  Node * root;
} AVLTree;



AVLTree * newAvlTree(int (*compareItems)(void *, void *), int (*compareIdItem)(int, void *)){
  AVLTree * new;

  new = malloc(sizeof(AVLTree));
  new->n_elements = 0;
  new->compareItems = compareItems;
  new->compareIdItem = compareIdItem;
  new->root = NULL;

  return new;
}

static int max(int a, int b){
  return a > b ? a : b;
}

static int getHeight(Node * node){
  return node == NULL ? 0 : node->height;
}

static Node * rightRotate(Node * y){
  Node * x = y->left;
  Node * T2 = x->right;

  // Perform rotation
  x->right = y;
  y->left = T2;

  // Update heights
  y->height = max(getHeight(y->left), getHeight(y->right))+1;
  x->height = max(getHeight(x->left), getHeight(x->right))+1;

  // Return new root
  return x;
}

static Node * leftRotate(Node * x){
  Node * y = x->right;
  Node * T2 = y->left;

  // Perform rotation
  y->left = x;
  x->right = T2;

  //  Update heights
  x->height = max(getHeight(x->left), getHeight(x->right))+1;
  y->height = max(getHeight(y->left), getHeight(y->right))+1;

  // Return new root
  return y;
}

static Node * insertNode(Node * root, Node * new_node, int (*compareItems)(void *, void *)){
  int balance;

  if(root == NULL){
    return new_node;
  }

  // Inserção de acordo com uma binary search tree (Nodes inferiores à esquerda e superiores à direita)
  if (compareItems(new_node->item, root->item) == SMALLER) {
    root->left = insertNode(root->left, new_node, compareItems);
  } else {
    root->right = insertNode(root->right, new_node, compareItems);
  }

  // Calcula informações para confirmar se a árvore está balanceada
  root->height = max(getHeight(root->left), getHeight(root->right)) + 1;
  balance = getHeight(root->left) - getHeight(root->right);

  // Faz as rotações necessárias para balancear a árvore
  if (balance > 1 && compareItems(new_node->item, root->left->item) == SMALLER){
    //Left Left
    return rightRotate(root);
  }
  if (balance < -1 && compareItems(new_node->item, root->right->item) == GREATER) {
    // Right Right
    return leftRotate(root);
  }
  if (balance > 1 && compareItems(new_node->item, root->left->item) == GREATER) {
    // Left Right
    root->left =  leftRotate(root->left);
    return rightRotate(root);
  }
  if (balance < -1 && compareItems(new_node->item, root->right->item) == SMALLER ) {
    // Right Left
    root->right = rightRotate(root->right);
    return leftRotate(root);
  }

  // Não é necessário alterações: Devolve o root original
  return root;
}

void insertAvlTree(AVLTree * tree, void * item){
  Node * new_node;

  new_node = malloc(sizeof(Node));
  new_node->item = item;
  new_node->height = 1;
  new_node->right = NULL;
  new_node->left = NULL;

  if (tree->n_elements == 0){
    tree->root = new_node;
  } else {
    tree->root = insertNode(tree->root, new_node, tree->compareItems);
  }

  tree->n_elements += 1;
}

void * findAvlTree(AVLTree * tree, int identifier){
  Node * aux;

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

static void freeNodes(Node * node, void (*freeItem)(void *)){
  if (node != NULL){
    freeNodes(node->left, freeItem);
    freeNodes(node->right, freeItem);
    freeItem(node->item);
    free(node);
  }
}

void freeAvlTree(AVLTree * tree, void (*freeItem)(void *)){
  if (tree != NULL){
    freeNodes(tree->root, freeItem);
    free(tree);
  }
}


void printAvlRecursive(Node * root, void (*printItem)(FILE *, int, int, void *), FILE * output_file, int x, int y){
  if (root != NULL){
    printAvlRecursive(root->left, printItem, output_file, x, y);
    printItem(output_file, x, y, root->item);
    printAvlRecursive(root->right, printItem, output_file, x, y);
  }
}

void printAvlTree(AVLTree * tree, void (*printItem)(FILE *, int, int, void *), FILE * output_file, int x, int y){
  if(tree != NULL){
    printAvlRecursive(tree->root, printItem, output_file, x, y);
  }
}

// *********************************************************************************************************************
// *********************************************************************************************************************








//                                            MAIN FILE
// *********************************************************************************************************************
// *********************************************************************************************************************
int SIZE;

#define cenas BADJORAS + 1

#define ALIVE 1
#define DEAD 0

#define TOP 0
#define BOTTOM 1
#define X_RIGHT 2
#define X_LEFT 3
#define Y_RIGHT 4
#define Y_LEFT 5

#define GET_INDEX(x, y) (SIZE)*(y) + (x)

#define checkVectorLimits(c) ((c) == SIZE ? 0 : ((c) == -1 ? SIZE-1 : (c)))

#define POV_Z(i, z) ((i) == TOP     ? z+1 : ((i) == BOTTOM ? z-1 : (z)))
#define POV_X(i, x) ((i) == X_RIGHT ? x+1 : ((i) == X_LEFT ? x-1 : (x)))
#define POV_Y(i, y) ((i) == Y_RIGHT ? y+1 : ((i) == Y_LEFT ? y-1 : (y)))

#define GET_X(i, x) checkVectorLimits(POV_X(i, x))
#define GET_Y(i, y) checkVectorLimits(POV_Y(i, y))
#define GET_Z(i, z) checkVectorLimits(POV_Z(i, z))

const int POV[] = {BOTTOM, TOP, X_LEFT, X_RIGHT, Y_LEFT, Y_RIGHT};

typedef struct cell{
  int x, y, z;
  int neighbors[6];
} Cell;


//Candidatas a macros!!
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


//------------------------------------------------------------------------------

void insertNewCellInAVL(AVLTree * tree, int x, int y, int z){
  Cell * new_cell;

  new_cell = malloc(sizeof(Cell));

  new_cell->x = x;
  new_cell->y = y;
  new_cell->z = z;
  memset(new_cell->neighbors, -1, 6*sizeof(int)); //Inicializa a memória com -1

  insertAvlTree(tree, new_cell);
}

void printCell(FILE * output_file, int x, int y, void * cell){
  if (cell != NULL){ //DEBUG
    fprintf(output_file, "%d %d %d\n", x, y, ((Cell *)cell)->z);
  }
}

AVLTree ** LoadMap(char * file_name){
  FILE * input_file;
  AVLTree ** map;
  int x, y, z;

  input_file = fopen(file_name, "r");

  if (input_file == NULL || fscanf(input_file, "%d", &SIZE) != 1){
    //printf("Erro no ficheiro de entrada\n");
    exit(0);
  }


  map = calloc(SIZE*SIZE, sizeof(AVLTree *)); //size^2 porque este vetor representa as coordenadas x e y simultaneamente (matriz)
  //map é logo inicializado a NULL porque é usado o calloc

  while(fscanf(input_file, "%d %d %d", &x, &y, &z) == 3){
    if (map[GET_INDEX(x,y)] == NULL){
      //Inicializa a AVL
      map[GET_INDEX(x,y)] = newAvlTree(compareItems, compareIdItem);
    }

    //O elemento a inserir na AVL é a coordenada Z
    insertNewCellInAVL(map[GET_INDEX(x,y)], x, y, z);
  }

  fclose(input_file);

  return map;
}



void addToNewGeneration(AVLTree ** next_map, int x, int y, int z, int nAliveNeighbors, int cell_state){
  if (nAliveNeighbors == 2 || nAliveNeighbors == 3 || (nAliveNeighbors == 4 && cell_state == ALIVE)){ // Condição para a proxima celula estar viva

    // Verifica se árvore corresponde às coordenadas x e y da celula já está criada
    if (next_map[GET_INDEX(x,y)] == NULL){
        next_map[GET_INDEX(x,y)] = newAvlTree(compareItems, compareIdItem); //Inicializa a AVL
    }

    // Aloca nova célula
    insertNewCellInAVL(next_map[GET_INDEX(x,y)], x, y, z);
  }
}

int visitNeighborsDead(AVLTree ** map, int x, int y, int z, int pov_caller){
  int i;
  int count = 1; //Uma célula morta só é chamada por uma célula viva!
  Cell *  neighbor_cell;

  for(i=0; i<6; i++){ //Vizita os 6 vizinhos de uma célula
    if (i != pov_caller){
      neighbor_cell = findAvlTree(map[GET_INDEX(GET_X(i, x), GET_Y(i, y))], GET_Z(i, z));

      if (neighbor_cell != NULL){
        // A célula vizinha está viva!!
        // Informa o vizinho sobre o estado da célula que o chamou
        neighbor_cell->neighbors[POV[i]] = DEAD;
        count ++; // A célula vizinha foi encontrada, por isso está viva
      }
    }
  }

  return count;
}



int visitNeighborsAlive(AVLTree ** map, AVLTree ** next_map, Cell * current_cell){
  int i, nAliveNeighbors;
  int count = 0;
  Cell *  neighbor_cell;

  for(i=0; i<6; i++){ //Vizita os 6 vizinhos de uma célula
    if (current_cell->neighbors[i] == -1){
      neighbor_cell = findAvlTree(map[GET_INDEX(GET_X(i, current_cell->x), GET_Y(i, current_cell->y))], GET_Z(i, current_cell->z));

      if (neighbor_cell != NULL){
        current_cell->neighbors[i] = ALIVE;

        // A célula vizinha está viva!!
        // Informa o vizinho sobre o estado da célula que o chamou
        neighbor_cell->neighbors[POV[i]] = ALIVE;
        count ++; // A célula vizinha foi encontrada, por isso está viva

      } else { //Avaliada a célula morta (Todas as celulas mortas que renascem têm pelo menos uma célula viva)
        current_cell->neighbors[i] = DEAD;

        nAliveNeighbors = visitNeighborsDead(map, current_cell->x, current_cell->y, current_cell->z, POV[i]);
        addToNewGeneration(next_map, current_cell->x, current_cell->y, current_cell->z, nAliveNeighbors, DEAD);
      }
    } else count += current_cell->neighbors[i];
  }

  return count;
}


void searchTree(AVLTree ** map, AVLTree ** next_map, Node * root){
  Cell * current_cell;
  int nAliveNeighbors;

  if (root != NULL){
    searchTree(map, next_map, root->left);
    searchTree(map, next_map, root->right);

    current_cell = root->item;
    nAliveNeighbors = visitNeighborsAlive(map, next_map, current_cell);
    addToNewGeneration(next_map, current_cell->x, current_cell->y, current_cell->z, nAliveNeighbors, ALIVE);
  }
}

void getNextGeneration(AVLTree ** map, AVLTree ** next_map){
  int xy;

  memset(next_map, 0, SIZE*SIZE*sizeof(AVLTree *)); //Inicializa o vetor next_map com o valor NULL

  for(xy= 0; xy<SIZE*SIZE; xy++){
    if (map[xy] != NULL){
      searchTree(map, next_map, map[xy]->root);
    }
  }
}

void freeMapAVLTrees(AVLTree ** map){
  int xy;

  for(xy=0; xy<SIZE*SIZE; xy++){
      freeAvlTree(map[xy], free);
  }
}

void writeOutput(char * file_name, AVLTree ** map){
  FILE * output_file;
  int x, y;

  output_file = fopen(file_name, "w");


  for(x=0; x<SIZE; x++){
    for(y=0; y<SIZE; y++){
      printAvlTree(map[GET_INDEX(x,y)], printCell, output_file, x, y);
    }
  }

  fclose(output_file);
}

int main(int argc, char *argv[]){
  int i, n_generations;
  AVLTree ** map, ** next_map, ** aux_pointer;
  char * file_name;

  if (argc != 3){
    //printf("Número argumentos não previsto\n");
    exit(0);
  }

  if (sscanf(argv[2],"%d", &n_generations) != 1){
    //printf("O segundo argumento não é um inteiro\n");
    exit(0);
  }

  file_name = malloc((strlen(argv[1])+2)*sizeof(char));
  strcpy(file_name, argv[1]);

  map = LoadMap(file_name);

  next_map = malloc(SIZE*SIZE*sizeof(AVLTree *));
  for (i=0; i<n_generations; i++){
    getNextGeneration(map, next_map);

    aux_pointer = map;
    map = next_map;
    next_map = aux_pointer;

    freeMapAVLTrees(next_map);
  }

  strcpy(file_name + strlen(file_name) - 2, "out\0"); //Substitui .in por .out


  writeOutput(file_name, map);

  freeMapAVLTrees(map);
  free(map);
  free(next_map);
  free(file_name);

  return 0;
}
