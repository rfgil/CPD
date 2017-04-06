#include <stdlib.h>
#include <stdio.h>

#include "avl_tree.h"

typedef struct node {
  void * item;
  int height;
  struct node * left;
  struct node * right;
} Node;

struct tree {
  int n_elements;
  int (*compareItems) (void *, void *);
  int (*compareIdItem) (int, void *);
  Node * root;
};

AVLTree * newAvlTree(int (*compareItems)(void *, void *), int (*compareIdItem)(int, void *)){
  AVLTree * new;

  new = malloc(sizeof(AVLTree));
  new->n_elements = 0;
  new->compareItems = compareItems;
  new->compareIdItem = compareIdItem;
  new->root = NULL;

  return new;
}

int countElementsAvlTree(AVLTree * tree){
  return tree->n_elements;
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

void * getAvlTreeRootNode(AVLTree * tree){
  return tree->root;
}

void * getAvlTreeLeftChildNode(void * node){
  return node == NULL ? NULL : ((Node *)node)->left;
}

void * getAvlTreeRightChildNode(void * node){
  return node == NULL ? NULL : ((Node *)node)->right;
}

void * getAvlTreeNodeItem(void * node){
  return node == NULL ? NULL : ((Node *)node)->item;
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
