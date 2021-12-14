
#include <stdio.h>
#include <stdlib.h>
#include "intLinkedList.h"

void appendI(intNodeI *target, int data){
	intNodeI *nu = (intNodeI*)malloc(sizeof(intNodeI));
	nu->data = data;
	nu->next = NULL;
	target->next = nu;
}

void print_listI(intNodeI *root){
	intNodeI *curr = root;
	while(curr != NULL){
		printf("%d", curr->data);
		printf(" ");
		curr = curr->next;
	}
	printf("\n");
}

void free_listI(intNodeI *root){
	intNodeI *node = root;
	while(node->next != NULL){
		intNodeI *temp = node;
		node = node->next;
		free(temp);
	}
	root = NULL;
}

void delete_nodeI(NodeI *it){
	free_listI(it->data);
	free(it);
}


NodeI* rmgetnodeI(intLinkedList *list, int i){
                if(i!=0){
			NodeI *curr = list->root;
                	int j;
                	NodeI *prev;
                	for(j=0; j<i; j++){
                        	prev = curr;
                        	curr = curr->next;
                	}
                	prev->next=curr->next;
                	NodeI *temp = curr;
			return temp;
		}
		else{
			NodeI *temp=list->root;
			list->root=temp->next;
			return temp;
		}
}


void List_printI(intLinkedList *list){
	NodeI *curr = list->root;
	int i = 0;
	while(curr != NULL){
		printf("Group %d:", i);
		print_listI(curr->data);
		curr = curr->next;
		i++;
	}
}

void insertendI(intLinkedList *list, NodeI *data){
	if(list->root!=NULL){
		NodeI *temp = list->root;
		while(temp-> next!=NULL){
			temp=temp->next;
		}
		NodeI *nu=(NodeI*)malloc(sizeof(NodeI));
		temp->next=nu;
		nu->data=data->data;
		nu->next=NULL;
	}
	else{
		NodeI *nu =(NodeI*)malloc(sizeof(NodeI));
		nu->data=data->data;
		list->root=nu;
		nu->next=NULL;
	}
}

void insertI(intLinkedList *list, int node_position, int data){
	int length = list->num_elems + 1;
	if(node_position > length){
		printf("Node insertion requested out of order.\n");
		return;
	}
	else{
 		NodeI *temp = list->root;
		int pos = 0;
		//traverse to the position
		while(pos < node_position){
			if(temp->next != NULL)
        			temp = temp->next;
      			else{
				temp->next = (NodeI*)malloc(sizeof(NodeI));
				temp->next->data = NULL;
				temp->next->next = NULL;
        			temp = temp->next;
      			}
      			pos++;
    		}
    		if(temp->data != NULL){
      			intNodeI *hop = temp->data;
      			while(hop->next != NULL)
        			hop = hop->next;
      			appendI(hop, data);
    		}
    		else{
      			intNodeI *nu = (intNodeI*)malloc(sizeof(intNodeI));
      			nu->data = data;
			nu->next = NULL;
      			//set Node's data to the newly created CharNode
      			temp->data = nu;
    		}
     		list->num_elems = length;
  	}
}

void free_Ilist(struct intLinkedList *list){
  	NodeI *temp = list->root;
  	while(temp->next != NULL){ // iterate through all Nodes in list
    		NodeI curr = *temp;         // Save node hopping along
    		temp = temp->next;
    		free_listI(curr.data);  // Free the list stored at the curr node
    		list->num_elems -= 1;
  	}
  	list = NULL;
}
