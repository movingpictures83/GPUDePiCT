#ifndef LINKEDLIST_H
#define LINKEDLIST_H

#include <stdio.h>

typedef struct CharNode{
	char            *data;
	struct CharNode *next;
}CharNode;

void append(CharNode *target, char *data);

void print_list(CharNode *root, int seqlength);

void free_listC(CharNode *root);

typedef struct Node{
	CharNode     *data;
	struct Node *next;
}Node;

void delete_node(Node *it);

typedef struct LinkedList{
	Node *root;
	int num_elems;
}LinkedList;

Node* rmgetnode(LinkedList *list, int i);


void List_print(LinkedList *list, int seqlength);

void swap(LinkedList *list, int to, int from);

void insertend(LinkedList *list, Node *data);

void insert(LinkedList *list, int node_position, char *data);

void free_list(struct LinkedList *list);

#endif
