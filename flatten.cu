#include "flatten.h"


void flatten(Node *curr, char *align, int *grouping, int rows, int seq_length){
	/* Flatten linked list and allocate it on the GPU if necessary*/
	int group_size =0, cur_group = 0, cur_position = 0;
	while(curr != NULL){
		CharNode *current = curr->data;
		while(current != NULL){
			int a;
			for (a = 0; a < seq_length; a++)
				align[cur_position*seq_length+a] = (current->data)[a];
			if (grouping!=NULL) group_size++;
			current = current->next;
			cur_position++;
		}
			if (grouping!=NULL){
				grouping[cur_group]=cur_position-1;
				cur_group +=1;
			}
			curr = curr->next;
	}
	if(grouping!=NULL){
		while(cur_group < rows){
			grouping[cur_group]=-1;
			cur_group++;
		}
	}
}
