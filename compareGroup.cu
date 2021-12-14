#include "compareGroup.h"

int compareGroup(int* blockScore, LinkedList *h_aligned, int *done, int* r, intLinkedList *tracer){
	int max_similarity=0, max_i, max_j, i;
	for(i=0; i<*r; i++){ // find max i & max j
		int similarity=0, itmp, jtmp=-1, j;
		for(j=i+1; j<*r; j++){
			if(*(blockScore+*r*i+j)>0){
				done[i]++;
				done[j]++;
			}
			if(*(blockScore+*r*i+j) > similarity){
				itmp=i;
				jtmp=j;
				similarity=*(blockScore+*r*i+j);
			}
		}
		/*if(similarity!=0){
			done[i]++;
			if(jtmp != -1){
				done[jtmp]++;
			}
		}*/
		if(similarity>max_similarity){
			max_i=itmp;
			max_j=jtmp;
			max_similarity=similarity;
			//printf("similarity: %i\n", similarity);
		}
//printf("%i\n", i);
	}
//printf("%d", max_similarity);
	if(max_similarity != 0)
		printf("Merge groups at indicies (%d, %d) with a score of: %d\n", max_i, max_j, max_similarity);
	else
		return -1;//break;

	Node *cur = h_aligned->root;
	NodeI *curT = tracer->root;
	Node *previous;
	NodeI *previousT;
	int cur_pos = 1;

	// ITERATE TO MAX J
	while(cur_pos < max_j+1){
		previous = cur;
		previousT = curT;
		cur = cur->next;
		curT = curT->next;
		cur_pos += 1;
	}

//printf("before charNode\n");
	CharNode *group = cur->data;
	intNodeI *groupT = curT->data;
	previous->next = cur->next;
	previousT->next = curT->next;
	free(cur);
	free(curT);
	h_aligned->num_elems--;
	tracer->num_elems--;
	*r = h_aligned->num_elems;

//printf("before iterate to max\n");
	// ITERATE TO MAX I & GROUP
	cur_pos = 1;
	Node *current = h_aligned->root;
	NodeI *currentT = tracer->root;
//printf("before while\n");
	while(cur_pos < max_i+1){
		current = current->next;
		currentT = currentT->next;
		cur_pos += 1;
	}
//printf("after while\n");
	CharNode *c = current->data;
	intNodeI *cT = currentT->data;
	while(c->next != NULL){
		c = c->next;
		cT = cT->next;
	}
	c->next = group;
	cT->next = groupT;
	return max_j;
}

