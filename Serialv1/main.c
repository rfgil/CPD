#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Estruturas/avl_tree.h"


#define ALIVE 1
#define DEAD 0

#define TOP 0
#define BOTTOM 1
#define X_RIGHT 2
#define X_LEFT 3
#define Y_RIGHT 4
#define Y_LEFT 5

typedef struct cell{
  int z;
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
#define GET_INDEX(x, y, size) (size)*(y) + (x)

#define checkVectorLimits(c, size) ((c) == size ? 0 : ((c) == -1 ? size-1 : (c)))

#define POV_Z(i, z, size) ((i) == TOP     ? z+1 : ((i) == BOTTOM ? z-1 : (z)))
#define POV_X(i, x, size) ((i) == X_RIGHT ? x+1 : ((i) == X_LEFT ? x-1 : (x)))
#define POV_Y(i, y, size) ((i) == Y_RIGHT ? y+1 : ((i) == Y_LEFT ? y-1 : (y)))

#define GET_X(i, x, size) checkVectorLimits(POV_X(i, x, size), size)
#define GET_Y(i, y, size) checkVectorLimits(POV_Y(i, y, size), size)
#define GET_Z(i, z, size) checkVectorLimits(POV_Z(i, z, size), size)


const int POV[] = {BOTTOM, TOP, X_LEFT, X_RIGHT, Y_LEFT, Y_RIGHT};

/*
void checkVectorLimits(int * coordinate, int size){
  if ((*coordinate) == size){
    (*coordinate) = 0;

  } else if((*coordinate) == -1){
    (*coordinate) = size-1;
  }
}
*/
/*
void getPOV(int neighbor_pos, int size, int x, int y, int z, int * new_x, int * new_y, int * new_z, int * pov){
  switch (neighbor_pos) {
    case TOP:
      (*new_z) = z + 1;
      (*new_z) = checkVectorLimits(new_z, size);
      (*pov) = BOTTOM;
      break;

    case BOTTOM:
      (*new_z) = z - 1;
      checkVectorLimits(new_z, size);
      (*pov) = TOP;
      break;

    case X_RIGHT:
      (*new_x) = x + 1;
      checkVectorLimits(new_x, size);
      (*pov) = X_LEFT;
      break;

    case X_LEFT:
      (*new_x) = x - 1;
      checkVectorLimits(new_x, size);
      (*pov) = X_RIGHT;
      break;

    case Y_RIGHT:
      (*new_y) = y + 1;
      checkVectorLimits(new_y, size);
      (*pov) = Y_LEFT;
      break;

    case Y_LEFT:
      (*new_y) = y - 1;
      checkVectorLimits(new_y, size);
      (*pov) = Y_RIGHT;
      break;

    default:
      return;
  }

  Se estou a visitar o vizinho de cima, eu próprio correspondo ao seu vizinho de baixo:
  0 - Top       -> 1 - Bottom
  1 - Bottom    -> 0 - Top
  2 - X Right   -> 3 - X Left
  3 - X Left    -> 2 - X Right
  4 - Y Right   -> 5 - Y Left
  5 - Y Left    -> 4 - Y Right

}

int getIndex(int x, int y, int size){
  return (size*y + x);
}
*/



//------------------------------------------------------------------------------

void printCell(FILE * output_file, int x, int y, void * cell){
  if (cell != NULL){ //DEBUG
    fprintf(output_file, "%d %d %d\n", x, y, ((Cell *)cell)->z);
  }
}

AVLTree ** LoadMap(char * file_name, int * size){
  FILE * input_file;
  Cell * new_cell;
  AVLTree ** map;
  int x, y, z;

  input_file = fopen(file_name, "r");

  if (input_file == NULL || fscanf(input_file, "%d", size) != 1){
    //printf("Erro no ficheiro de entrada\n");
    exit(0);
  }


  map = calloc((*size)*(*size), sizeof(AVLTree *)); //size^2 porque este vetor representa as coordenadas x e y simultaneamente (matriz)
  //map é logo inicializado a NULL porque é usado o calloc

  /*
  for(i = 0; i<(*size)*(*size); i++){
    map[i] = NULL;
  }
  */

  while(fscanf(input_file, "%d %d %d", &x, &y, &z) == 3){
    if (map[GET_INDEX(x,y,*size)] == NULL){
      //Inicializa a AVL
      map[GET_INDEX(x,y,*size)] = newAvlTree(compareItems, compareIdItem);
    }

    //O elemento a inserir na AVL é a coordenada Z
    new_cell = malloc(sizeof(Cell));
    new_cell->z = z;
    memset(new_cell->neighbors, -1, 6*sizeof(int)); //Inicializa a memória com -1

    insertAvlTree(map[GET_INDEX(x,y,*size)], new_cell);
  }

  fclose(input_file);

  return map;
}


int visitNeighbors(AVLTree ** map, int size, int * neighbors, int x, int y, int z, int is_alive){
  int i;
  int count = 0;
  Cell *  n_cell;

  for(i=0; i<6; i++){ //Vizita os 6 vizinhos de uma célula
    if (neighbors[i] == -1){
      n_cell = findAvlTree(map[GET_INDEX(GET_X(i, x, size), GET_Y(i, y, size), size)], GET_Z(i, z, size));

      if (n_cell != NULL){
        n_cell->neighbors[POV[i]] = is_alive; // Informa o vizinho sobre o estado da célula que o chamou
        count ++; // A célula vizinha foi encontrada, por isso está viva
      }
    } else count += neighbors[i];
  }

  return count;
}

void getNextGeneration(AVLTree ** map, AVLTree ** next_map, int size){
  int x, y, z;

  int neighbors[6];
  int count;

  Cell * new_cell, * current_cell;


  memset(next_map, 0, size*size*sizeof(AVLTree *)); //Inicializa o vetor next_map com o valor NULL

  for(xy= 0; xy<size*size; xy++){
    
  }


  for(x=0; x<size; x++){
    for(y=0; y<size; y++){
      for(z=0; z<size; z++){

        current_cell = (Cell *) findAvlTree(map[GET_INDEX(x, y, size)], z);

        if(current_cell != NULL){
          count = visitNeighbors(map, size, current_cell->neighbors, x, y, z, ALIVE);
        } else {
          memset(neighbors, -1, 6*sizeof(int)); // Inicializa neighbors com -1
          count = visitNeighbors(map, size, neighbors, x, y, z, DEAD);
        }

        if (count == 2 || count == 3 || (count == 4 && current_cell != NULL)){ // Condição para a proxima celula estar viva
          if (next_map[GET_INDEX(x,y,size)] == NULL){
            next_map[GET_INDEX(x,y,size)] = newAvlTree(compareItems, compareIdItem); //Inicializa a AVL
          }

          new_cell = malloc(sizeof(Cell));
          new_cell->z = z;
          memset(new_cell->neighbors, -1, 6*sizeof(int)); //Inicializa a memória com -1

          insertAvlTree(next_map[GET_INDEX(x,y,size)], new_cell);
        }
      }
    }
  }
}



void freeMapAVLTrees(AVLTree ** map, int size){
  int xy;

  for(xy=0; xy<size*size; xy++){
      freeAvlTree(map[xy], free);
  }
}



void writeOutput(char * file_name, AVLTree ** map, int size){
  FILE * output_file;
  int x, y;

  output_file = fopen(file_name, "w");


  for(x=0; x<size; x++){
    for(y=0; y<size; y++){
      printAvlTree(map[GET_INDEX(x,y,size)], printCell, output_file, x, y);
    }
  }

  fclose(output_file);
}

int main(int argc, char *argv[]){
  int size, i, n_generations;
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

  map = LoadMap(file_name, &size);

  next_map = malloc(size*size*sizeof(AVLTree *));
  for (i=0; i<n_generations; i++){
    getNextGeneration(map, next_map, size);

    aux_pointer = map;
    map = next_map;
    next_map = aux_pointer;

    freeMapAVLTrees(next_map, size);
  }

  strcpy(file_name + strlen(file_name) - 2, "out\0"); //Substitui .in por .out


  writeOutput(file_name, map, size);

  freeMapAVLTrees(map, size);
  free(map);
  free(next_map);
  free(file_name);

  return 0;
}
