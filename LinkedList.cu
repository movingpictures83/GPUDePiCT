
#include <stdio.h>
#include <stdlib.h>
#include "LinkedList.h"

void append(CharNode *target, char *data){
	CharNode *nu = (CharNode*)malloc(sizeof(CharNode));
	nu->data = data;
	nu->next = NULL;
	target->next = nu;
}

void print_list(CharNode *root, int seqlength){
	CharNode *curr = root; int i;
	while(curr != NULL){
		for(i = 0; i < seqlength; i++)
			printf("%c", (curr->data)[i]);
		printf("  ");
		curr = curr->next;
	}
	printf("\n");
}

void free_listC(CharNode *root){
	CharNode *node = root;
	while(node->next != NULL){
		CharNode *temp = node;
		node = node->next;
		free(temp->data);
	}
	root = NULL;
}

void delete_node(Node *it){
	free_listC(it->data);
	free(it);
}


Node* rmgetnode(LinkedList *list, int i){
                if(i!=0){
			Node *curr = list->root;
                	int j;
                	Node *prev;
                	for(j=0; j<i; j++){
                        	prev = curr;
                        	curr = curr->next;
                	}
                	prev->next=curr->next;
                	Node *temp = curr;
			return temp;
		}
		else{
			Node *temp=list->root;
			list->root=temp->next;
			return temp;
		}
}


void List_print(LinkedList *list, int seqlength){
	Node *curr = list->root;
	int i = 0;
	while(curr != NULL){
		printf("NEXT %d:\n", i);
		print_list(curr->data, seqlength);
		curr = curr->next;
		i++;
	}
}

void swap(LinkedList *list, int to, int from){
	Node *from_n = list->root;
	Node *to_n   = list->root;
	CharNode *temp;  CharNode *tmp;
	int cur_pos = 1;
	while(cur_pos < from){
 		from_n = from_n->next;
		cur_pos++;
	}
	temp = from_n->data;
	cur_pos = 1;
	while(cur_pos < to){
 		to_n = to_n->next;
		cur_pos++;
	}
	tmp = to_n->data;
	/* make the switch */
	from_n->data =    tmp;
	to_n->data   =   temp;
}

void insertend(LinkedList *list, Node *data){
	if(list->root!=NULL){
		Node *temp = list->root;
		//while(temp->next!=NULL){
		while(temp!=NULL){
			temp=temp->next;
		}
		Node *nu=(Node*)malloc(sizeof(Node));
		temp->next=nu;
		nu->data=data->data;
		nu->next=NULL;
	}
	else{
		Node *nu =(Node*)malloc(sizeof(Node));
		nu->data=data->data;
		list->root=nu;
		nu->next=NULL;
	}
}

void insert(LinkedList *list, int node_position, char *data){
	int length = list->num_elems + 1;
	if(node_position > length){
		printf("Node insertion requested out of order.\n");
		return;
	}
	else{
 		Node *temp = list->root;
		int pos = 1;
		//traverse to the position
		while(pos < node_position){
			if(temp->next != NULL)
        			temp = temp->next;
      			else{
				temp->next = (Node*)malloc(sizeof(Node));
				temp->next->data = NULL;
				temp->next->next = NULL;
        			temp = temp->next;
      			}
      			pos++;
    		}
    		if(temp->data != NULL){
      			CharNode *hop = temp->data;
      			while(hop->next != NULL) {
        			hop = hop->next;
			}
      			append(hop, data);
    		}
    		else{
      			CharNode *nu = (CharNode*)malloc(sizeof(CharNode));
      			nu->data = data;
			nu->next = NULL;
      			//set Node's data to the newly created CharNode
      			temp->data = nu;
    		}
     		list->num_elems = length;
  	}
}

void free_list(struct LinkedList *list){
  	Node *temp = list->root;
  	while(temp->next != NULL){ // iterate through all Nodes in list
    		Node curr = *temp;         // Save node hopping along
    		temp = temp->next;
    		free_listC(curr.data);  // Free the list stored at the curr node
    		list->num_elems -= 1;
  	}
  	list = NULL;
}

