
#include "compareGroupK.h"

int compareGroupK(int* blockScore, LinkedList *h_aligned, int* done, int* r, int* moreCents, int* cents, int* numCents, intLinkedList *tracer){
	//flag to indicate first run through while
	int flag=1;
	//runs until all distinct centroids are found
	while(1){
		int i, iflag=0;
		int max_similarity=0;
		for(i=0; i<*r; i++){ // find max i & max j
			//increments i to a spot previous centroids have a score of 0 at this point
			while(*(moreCents+i)!=-1){
				i++;
				//at end break out of while and for i
				if(i==*r){
					iflag=1;
					break;
				}
			}
			if(iflag==1){
				break;
			}
			int similarity=0, itmp, jtmp=-1, j, jflag=0;
			for(j=i+1; j<*r; j++){
				//increments j to a spot previous centroids have a score of 0 at this point
				while(*(moreCents+j)!=-1) {
					j++;
					//at end break out of while and for j
					if(j==*r){
						jflag=1;
						break;
					}
				}
				if(jflag==1)
					break;
				if(*(blockScore+*r*i+j) > 0&&flag==1){
					(*(done+i))++;
					(*(done+j))++;
				}
				//updates similarity to highest blockScore[i][?]
				if(*(blockScore+*r*i+j) > similarity){
					itmp=i;
					jtmp=j;
					similarity=*ptrMath2D(blockScore, i, j, *r);
				}
			}
			//no similar groups to i and on first iteration
			//done[i] and done[jtmp] are incremented
			/*if(similarity!=0&&flag==1){
				(*(done+i))++;
				if(jtmp != -1){
					(*(done+jtmp))++;
				}
			}*/
			//updates max_similarity, puts values in cents
			if(similarity>max_similarity){
				//max_i=itmp;
				//max_j=jtmp;
				*(cents+*numCents*2+0)=itmp;
				*(cents+*numCents*2+1)=jtmp;
				max_similarity=similarity;
			}
		}
		//increments moreCents where the centroid has a score of nonzero and the centroid itself
		if(max_similarity != 0){
			flag=0;
			printf("Merge groups at indicies (%d, %d) with a score of: %d\n", *(cents+*numCents*2+0), *(cents+*numCents*2+1), max_similarity);
			for(i=0; i<*r; i++){
				if(*(blockScore+*r*(*(cents+*numCents*2+0))+i)>0||*(blockScore+*r*i+(*(cents+*numCents*2+0)))>0){
					(*(moreCents+i))++;
				}
				else if(i==*(cents+*numCents*2+0)){
					(*(moreCents+i))++;
				}
			}
		}
		//no Centroids, max_similarity is 0, done with main while
		else if(*numCents==0)
			return -1;//break;
		//max_similarity is 0, done finding centroids
		else
			return 1;
		//merges centroid and match
		Node *cur = h_aligned->root;
		NodeI *curT = tracer->root;
		Node *previous;
		NodeI *previousT;
		int cur_pos = 1;
		// ITERATE TO MATCH
		while(cur_pos < *(cents+*numCents*2+1)+1){
			//accounts for centroids already removed in list
			for(i=0; i<*numCents; i++){
				//if(cur_pos==*(cents+i*2+0)||cur_pos==*(cents+i*2+1)){
				if(cur_pos==*(cents+i*2+1)){
					cur_pos++;
				}
			}
			previous = cur;
			previousT = curT;
			cur = cur->next;
			curT = curT->next;
			cur_pos += 1;
		}
		//removes match
		CharNode *group = cur->data;
		intNodeI *groupT = curT->data;
		//previous->next = cur->next;
		//previousT->next = curT->next;
			if(cur == h_aligned->root){
				h_aligned->root = cur->next;
				tracer->root = curT->next;
			}
			else{
                		previous->next = cur->next;
				previousT->next = curT->next;
			}
		free(cur);
		free(curT);
		h_aligned->num_elems--;
		tracer->num_elems--;
		//(*r) = h_aligned->num_elems;

		// ITERATE TO CENTROID & GROUP
		cur_pos = 1;
		Node *current = h_aligned->root;
		NodeI *currentT = tracer->root;
		while(cur_pos < *(cents+*numCents*2+0)){
			//accounts for centroids already removed in list
			for(i=0; i<*numCents; i++){
				//if(cur_pos==*(cents+i*2+0)||cur_pos==*(cents+i*2+1)){
				if(cur_pos==*(cents+i*2+1)) {
					cur_pos++;
				}
			}
			current = current->next;
			currentT = currentT->next;
			cur_pos += 1;
		}
		//groups centroid and match
		CharNode *c = current->data;
		intNodeI *cT = currentT->data;;
		while(c->next != NULL){
			c = c->next;
			cT = cT->next;
		}
		c->next = group;
		cT->next = groupT;
		//increments number of centroids
		(*numCents)++;
	}
}

void fuzzyK(int* blockScore, LinkedList *h_aligned, int* r, int* numCents, int* fuzzymatch, int* cents, int* fuzziness, int* fuzzies, intLinkedList *tracer, int* tolerance){
	int m, f;
			/*printf("**********************************************************************************\n");
			printf("* OUR LIST BEFORE ANYTHING                                                        \n");
			List_print(h_aligned, 1126);
			printf("**********************************************************************************\n");*/
	//printf("running fuzzy\n");
	//each iteration populates fuzzymatch[centroid][m]
	for(m=0; m<*fuzziness; m++){
		//printf("first for\n");
		int k;
		//finds a fuzzy match for each centroid
		for(k=0; k< *numCents; k++){
			//printf("second for k:%d numCents:%d\n", k, *numCents);
			int i;
			int max_simlow=0, match=-1, cent, breakflag;
			//finds the highest degree of membership
			for(i=0; i< *numCents; i++){
				breakflag=0;
				//printf("numcent for i:%d numCents:%d\n", i, *numCents);
				//checks that no fuzzymatch has been found for this spot
				while(*ptrMath2D(fuzzymatch, i, m, *fuzziness)!=-1){
					i++;
					//printf("i increased\n");
					if(i==*numCents){
						breakflag=1;
						break;
					}
				}
				if(breakflag==1)
					break;
				int j;
				int simlow=0, bposi, bposj, bposfuzzy, matchtmp, simlowtmp=0, fuzzyflag=0;
				//computes finds the highest low blockScore and the match for this centroid
				for(j=0; j<*r; j++){
					fuzzyflag=0;
					//printf("j:%d r:%d cent1:%d\n", j, *r, *ptrMath2D(cents, i, 0, 2));
					//assigns correct spot
					int b=0, bflag=0;
					while(b==0){
						b=1;
						int n;
						for(n=0; n<*numCents; n++){
                                        		//accounts for centroid matches already removed
                                        		if(j==*ptrMath2D(cents, n, 1, 2)||j==*ptrMath2D(cents, n, 0, 2)){
                                                		j++;
								b=0;
								//printf("j increased\n");
							}
                                        		//accounts for fuzzy matches already removed
                                        		for(f=0; f<m+1; f++){
								if(*ptrMath2D(fuzzymatch, n, f, *fuzziness)==-1)
									break;
                                                		else if(j==*ptrMath2D(fuzzymatch, n, f, *fuzziness)){
                                                        		j++;
									b=0;
									//printf("j increased\n");
								}
                                        		}
                                        		if(j >= *r){
								b=1;
								bflag=1;
								break;
							}
                                		}
					}
					if(bflag==1){
						//printf("breaking at bflag\n");
						break;
					}
					if(*ptrMath2D(cents, i, 0, 2)>j){
						bposi=*r*j+*ptrMath2D(cents, i, 0, 2);
						//printf("bposi:[%d][%d]\n", j, *ptrMath2D(cents, i, 0, 2));
					}
					else if(*ptrMath2D(cents, i, 0, 2)<j){
						bposi=*r*(*ptrMath2D(cents, i, 0, 2))+j;
						//printf("bposi:[%d][%d]\n", *ptrMath2D(cents, i, 0, 2), j);
					}
					//check that it's not 0, otherwise sets the flag (all matches must be non-zero)
					if(*(blockScore+bposi)!=0){
						simlowtmp=*(blockScore+bposi);
						//printf("simlowtmp1:%d\n", simlowtmp);
					}
					else
						fuzzyflag=1;
					if(*ptrMath2D(cents, i, 1, 2)>j){
                                	        bposj=*r*j+*ptrMath2D(cents, i, 1, 2);
						//printf("bposj:[%d][%d]\n", j, *ptrMath2D(cents, i, 1, 2));
					}
                                	else if(*ptrMath2D(cents, i, 1, 2)<j){
                                	        bposj=*r*(*ptrMath2D(cents, i, 1, 2))+j;
						//printf("bposj:[%d][%d]\n", *ptrMath2D(cents, i, 1, 2), j);
					}
					if(*(blockScore+bposj)!=0){
						//printf("nonzero\n");
						if(*(blockScore+bposj)<simlowtmp){
							simlowtmp=*(blockScore+bposj);
							//printf("simlowtmp2:%d\n", simlowtmp);
						}
					}
					else
						fuzzyflag=1;

					int f;
					//does the same as above for previously found fuzzy nodes
					for(f=1; f<m+1; f++){
						//printf("fuzzy for %d\n", f);
						if(*ptrMath2D(fuzzymatch, i, f-1, *fuzziness)!=-1){
							if(*ptrMath2D(fuzzymatch, i, f-1, *fuzziness)>j){
								bposfuzzy=*r*j+*ptrMath2D(fuzzymatch, i, f-1, *fuzziness);
								//printf("bposfuzzy:(%d, %d) i:%d\n", j, *ptrMath2D(fuzzymatch, i, f-1, *fuzziness), i);
							}
							else if(*ptrMath2D(fuzzymatch, i, f-1, *fuzziness)<j){
								bposfuzzy=(*r)*(*ptrMath2D(fuzzymatch, i, f-1, *fuzziness))+j;
								//printf("bposfuzzy:(%d, %d) i:%d\n", *ptrMath2D(fuzzymatch, i, f-1, *fuzziness), j, i);
							}
							//printf("in middle\n");
							if(*(blockScore+bposfuzzy)!=0 ) {
								//printf("Score:%d simlow: simlowtmp: %d\n", *(blockScore+bposfuzzy), simlowtmp);
								//printf("error\n");
								if(*(blockScore+bposfuzzy)<simlowtmp)
									simlowtmp=*(blockScore+bposfuzzy);
							}
							else
								fuzzyflag=1;
						}
						else
							break;
					}
					//printf("after fuzzy for\n");
					//updates sim low if all are non-zero scores
					if(fuzzyflag==0&&simlowtmp>simlow){
						//printf("simlow:%d\n", simlowtmp);
						simlow=simlowtmp;
						matchtmp=j;
						//printf("matchtmp:%d\n", matchtmp);
					}
				}



				float percentage = ((float)simlow/(float)(*ptrMath2D(blockScore, *ptrMath2D(cents, i, 0, 2), *ptrMath2D(cents, i, 1, 2), *r)));
				//updates max_simlow and match(highest match found)
				if(simlow>max_simlow && *tolerance<=100*percentage){
					cent=i;
					match=matchtmp;
					max_simlow=simlow;
					//printf("found! max_simlow:%d match:%d\n", max_simlow, match);
				}
			}
			//match found! updates fuzzymatch[cent][m] and increments fuzzies
			if(match!=-1){
				printf("Merging %d with centroid#: %d at position %d with low score of %d\n", match, cent, *ptrMath2D(cents, cent, 0, 2), max_simlow);
				*ptrMath2D(fuzzymatch, cent, m, *fuzziness)=match;
				(*fuzzies)++;
			}
			//no match found at highest degree of fuzziness desired, done
			else if(m==*fuzziness-1){
				//printf("return\n");
				return;
			}
			//no match found at this degree of fuzziness, may be some at higher level
			else{
				//printf("break\n");
				break;
			}
			//merges centroid and match
                	Node *cur = h_aligned->root;
			NodeI *curT = tracer->root;
                	Node *previous;
			NodeI *previousT;
                	int cur_pos = 0;
                	// ITERATE TO MATCH
			int monkeyflag=0;
			//printf("going to: %d\n", *ptrMath2D(fuzzymatch, cent, m, *fuzziness));
                	while(cur_pos < *ptrMath2D(fuzzymatch, cent, m, *fuzziness)){
				//accounts for centroids already removed in list
                	        //printf("current cur_pos:%d\n", cur_pos);
				int b=0;
				while(b==0){
					b=1;
					for(i=0; i<*numCents; i++){
						//accounts for centroid matches already removed
                	                	//printf("i:%d numCents:%d centmatch:%d\n", i, *numCents, *ptrMath2D(cents, i, 1, 2));
						if(cur_pos==*ptrMath2D(cents, i, 1, 2)){
							//printf("CURPOS IS %d, SO INCREMENTING: %d\n", cur_pos, *ptrMath2D(cents, i, 1, 2));
							//printf("pos++\n");
							b=0;
                	                        	cur_pos++;
							if(cur_pos >= *ptrMath2D(fuzzymatch, cent, m, *fuzziness)){
								b=1;
								monkeyflag=1;
								//printf("break1\n");
								break;
							}
						}
						//accounts for fuzzy matches already removed
						for(f=0; f<m+1; f++){
							if(cur_pos==*ptrMath2D(fuzzymatch, i, f, *fuzziness)&&(i!=cent || f!=m)){
							//printf("CURPOS IS %d, SO INCREMENTING: %d\n", cur_pos, *ptrMath2D(fuzzymatch, i, f, *fuzziness));
								//printf("pos++");
								b=0;
								cur_pos++;
								if(cur_pos >= *ptrMath2D(fuzzymatch, cent, m, *fuzziness)){
									b=1;
									monkeyflag=1;
									//printf("break2\n");
									break;
								}
							}
							else if(*ptrMath2D(fuzzymatch, i, f, *fuzziness)==-1&&(i!=cent || f!=m))
								break;
							else {
							//printf("CURPOS IS %d, NOT INCREMENTING: %d.  I: %d  CENT: %d, F: %d, M: %d\n", cur_pos, *ptrMath2D(fuzzymatch, i, f, *fuzziness), i, cent, f, m);

							}
						}
						//if(cur_pos >= *ptrMath2D(fuzzymatch, cent, m, *fuzziness)-1){
						//	printf("ever running?\n");
						//	b=1;
						//	monkeyflag=1;
						//	break;
						//}
                	        	}
				}
				if(monkeyflag==1)
					break;
                	        previous = cur;
				previousT = curT;
				//printf("JUMPING CUR\n");
                	        cur = cur->next;
				curT = curT->next;
                	        cur_pos += 1;
                	}
			/*printf("end of loop\n");
			printf("**********************************************************************************\n");
			printf("* OUR LIST AFTER FIRST MERGES                                                     \n");
			List_print(h_aligned, 1126);
			printf("**********************************************************************************\n");*/
			//printf("cur_pos: %d\n", cur_pos);
                	//removes match
			//printf("%d\n", (cur == NULL));
			//printf("%d\n", (cur->data == NULL));
                	CharNode *group = cur->data;
			//printf("group assigned\n");
			intNodeI *groupT = curT->data;
                	//previous->next = cur->next;
			//previousT->next = curT->next;
			if(cur == h_aligned->root){
				h_aligned->root = cur->next;
				tracer->root = curT->next;
			}
			else{
                		previous->next = cur->next;
				previousT->next = curT->next;
			}
                	free(cur);
			free(curT);
                	h_aligned->num_elems--;
			tracer->num_elems--;

                	// ITERATE TO CENTROID & GROUP
                	cur_pos = 0;
                	Node *current = h_aligned->root;
			NodeI *currentT = tracer->root;

			//printf("end of assignments\n");
			//printf("CENTROID GOING TO: %d", *ptrMath2D(cents, cent, 0, 2));
			monkeyflag=0;
			while(cur_pos < *ptrMath2D(cents, cent, 0, 2)){
				//printf("centroid while\n");
				int b=0;
				while(b==0){
					b=1;
					for(i=0; i<*numCents; i++){
			//printf("CENTROID %d SECOND VALUE: %d\n", i, *ptrMath2D(cents, i,  1, 2));
                                        	if(cur_pos==*ptrMath2D(cents, i, 1, 2)){
							//printf("CURPOS IS %d, SO INCREMENTING", cur_pos, *ptrMath2D(cents, i, 1, 2));
                                                	b=0;
							cur_pos++;
						}
                                        	for(f=0; f<m+1; f++){
							//printf("FUZZY VALUE (%d, %d): %d\n", i, f, *ptrMath2D(fuzzymatch, i, f, *fuzziness));
                                                	if(cur_pos==*ptrMath2D(fuzzymatch, i, f, *fuzziness)){
							//printf("CURPOS IS %d, FUZZY SO INCREMENTING", cur_pos, *ptrMath2D(fuzzymatch, i, f, *fuzziness));
								b=0;
                                                        	cur_pos++;
							}
                                        	}
                                        	if(cur_pos >= *ptrMath2D(cents, cent, 0, 2)){
							//printf("CURPOS IS %d, WHICH IS >= %d, SO BREAKING", cur_pos, *ptrMath2D(cents, cent, 0, 2));
                                                	b=1;
							monkeyflag=1;
							break;
						}
                                	}
				}
				if(monkeyflag!=0)
					break;
				//printf("JUMPING CURRENT\n");
                	        current = current->next;
				//print_list(current->data, 1126);
				currentT = currentT->next;
                	        cur_pos++;
				//printf("CUR POS END OF LOOP: %d\n", cur_pos);
                	}
			//printf("cur_pos: %d\n", cur_pos);
                	//groups centroid and match
                	CharNode *c = current->data;
			intNodeI *cT = currentT->data;
                	while(c->next != NULL){
                	        c = c->next;
				cT = cT->next;
			}
                	c->next = group;
			cT->next = groupT;
			/*printf("**********************************************************************************\n");
			printf("* OUR LIST CENTROID MERGE                                                         \n");
			List_print(h_aligned, 1126);
			printf("**********************************************************************************\n");*/
		}
	}
}
