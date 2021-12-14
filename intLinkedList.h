#ifndef INTLINKEDLIST_H
#define INTLINKEDLIST_H

#include <stdio.h>

typedef struct intNodeI{
	int             data;
	struct intNodeI *next;
}intNodeI;

void appendI(intNodeI *target, int data);

void print_listI(intNodeI *root);

void free_listI(intNodeI *root);

typedef struct NodeI{
	intNodeI     *data;
	struct NodeI *next;
}NodeI;

void delete_nodeI(NodeI *it);

typedef struct intLinkedList{
	NodeI *root;
	int num_elems;
}intLinkedList;

NodeI* rmgetnodeI(intLinkedList *list, int i);

void List_printI(intLinkedList *list);

void swapI(intLinkedList *list, int to, int from);

void insertendI(intLinkedList *list, NodeI *data);

void insertI(intLinkedList *list, int node_position, int data);

void free_Ilist(struct intLinkedList *list);
#endif
