#ifndef AVLTREE_H_INCLUDED
#define AVLTREE_H_INCLUDED

typedef struct tree AVLTree;

#define SMALLER 1
#define GREATER 2
#define EQUAL 3

AVLTree * newAvlTree(int (*compare)(void *, void *), int (*compareIdItem)(int, void *));
void insertAvlTree(AVLTree * tree, void * item);
void * findAvlTree(AVLTree * tree, int identifier);

void printAvlTree(AVLTree * tree, void (*printItem)(FILE *, int, int, void *), FILE * output_file, int x, int y);

/*
int countElementsAvlTree(AVLTree * tree);

void * getAvlTreeRootNode(AVLTree * tree);
void * getAvlTreeLeftChildNode(void * node);
void * getAvlTreeRightChildNode(void * node);
void * getAvlTreeNodeItem(void * node);
*/

void freeAvlTree(AVLTree * tree, void (*freeItem)(void *));

#endif
