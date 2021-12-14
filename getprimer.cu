#include "getprimer.h"

void getprimer(LinkedList* list, int *group, char *align, int groupsb4, char *simacids[NUMACIDS], char *codons[NUMACIDS], int seq_length, int min_primer_length){
	int index, lowerJ, n, i, j, k;
	for(i=0; i<list->num_elems; i++){
		//printf("=======================================================================================\n");
                //printf("PRIMERS FOR GROUP %d:\n", i+groupsb4);
                // Create a primer for the group

                if(i == 0)
                        lowerJ = 0;
                else
                        lowerJ = group[i-1]+1;
                //printf("Index range: %d to %d", lowerJ, group[i]);
                int start=0, count=0;
		int codonids[seq_length];
		for(k=0; k <seq_length; k++)
			codonids[k]=-1;
		int numprim = 0;
		for(k=0; k <seq_length; k++){
                        // Generate a string of acids
                        char acids[7];
                        acids[0] = '\0';
                        int size = 1;
                        for(j=lowerJ; j<= group[i]; j++){
                                //printf("K: %d   GROUP: %d\n", k, i+groupsb4);
                                char singleacid = align[j*seq_length+k];
				//printf("%c", singleacid);
                                // shift the acid into 'acids' string
                                int a = 0;
                                while(acids[a] < singleacid && acids[a] != '\0')
                                        a++;
                                // a is the first acid that is greater or equal to single acid
                                // if it is equal, nothing to do...
                                // if it is greater, shift!
                                if(acids[a] != singleacid){
                                        //printf("ELSE CASE A: %d \n", a);
                                        if(size != 7){
                                                n = size;
                                                //printf("BEFORE ALOOP %d %d\n", n, a);
                                                while(n !=  a){
                                                        acids[n] = acids[--n];
                                                }
                                                //printf("AFTER ALOOP\n");
                                                acids[a] = singleacid;
                                                //printf("ACIDS: %s", acids);
                                                size++;
					 }
                                }
                        }
                        int low=0, up=NUMACIDS-1, mid;
                        //printf("Binary search. ACIDS: %s\n", acids);
                        while(1){
                                mid = (up + low)/2;
                                //printf("MID IS: %d\n", mid);
                                //printf("Comparing %s and %s\n", simacids[mid], acids);
                                if( strcmp(simacids[mid], acids) == 0 ){
                                        //printf("Found: %s\n", acids);
					index = mid;
                                        break;
                                }
				// acid not found
                                else if(up < low){
					//printf("Not found: %s\n", acids);
                                        index = -1;
                                        break;
                                }
                                else if(strcmp(simacids[mid], acids) < 0)
                                        low = mid + 1;
                                else
                                        up = mid - 1;
                        }
                        //printf("Post binary search.\n");
			if(index != -1){
                                codonids[k]=index;
				count++;
				numprim++;
			}
                        else{
				if(count < min_primer_length){
					for(j=0; j < count; j++)
						codonids[j+start]=-1;
					numprim -= count;
				}
				start=k+1;
				count=0;
                        }
                } // end of k
		if(count < min_primer_length){
			for(j=0; j<count; j++)
				codonids[j+start]=-1;
			numprim -= count;
		}

		int endseq=0;
		int primercount=0;
		if (numprim != 0) {
		printf("=======================================================================================\n");
                printf("PRIMERS FOR GROUP %d:\n", i+groupsb4);
		//int endseq=0;
		//int primercount=0;
		for(k=0; k < seq_length; k++){
			if(codonids[k] != -1){
				if(endseq==0){
					printf("Primer %i @ %i:", primercount,k);
					primercount++;
					endseq=1;
				}
				printf("%s", codons[codonids[k]]);
			}
			else{
				if(endseq==1){
					printf("\n\n");
					endseq=0;
				}
			}
		}
                printf("\n");
		}
        } //end of primer creation
}

