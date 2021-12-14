#include "getNucPrimer.h"

void getNucPrimer(LinkedList* list, int *group, char *align, int groupsb4, char *simNucs[NUMNUCS], char *nucCodons[NUMNUCS], int seq_length, int min_primer_length){
	int index, lowerJ, n, i, j, k;
	for(i=0; i<list->num_elems; i++){
		//printf("=======================================================================================\n");
                //printf("PRIMERS FOR GROUP %d:\n", i+groupsb4);
                // Create a primer for the group
                if(i == 0)
                        lowerJ = 0;
                else
                        lowerJ = group[i-1]+1;
                //printf("Index range: %d to %d \n", lowerJ, group[i]);
                int start=0, count=0;
		int codonids[seq_length];
		for(k=0; k <seq_length; k++)
			codonids[k]=-1;
		int numprim=0;
		int count2 = 0;
		for(k=0; k <seq_length; k++){
                        // Generate a string of nucleotides
                        char nucs[5];
                        nucs[0] = '\0';
                        int size = 1;
                        for(j=lowerJ; j<= group[i]; j++){
                                //printf("K: %d GROUP: %d\n", k, i+groupsb4);
                                char singlenuc = align[j*seq_length+k];
                                //printf("SINGLENUC: %c", singlenuc);
                                // shift the nucleotide into 'nucs' string
                                int a = 0;
                                while(nucs[a] < singlenuc && nucs[a] >= 65)
                                        a++;
                                // a is the first nucleotide that is greater or equal to single nucleotide
                                // if it is equal, nothing to do...
                                // if it is greater, shift!
                                if(nucs[a] != singlenuc){
                                        //printf("ELSE CASE A: %d \n", a);
                                        if(size != 5){
                                                n = size;
                                                //printf("BEFORE ALOOP %d %d\n", n, a);
                                                while(n != a){
							nucs[n] = nucs[n-1];
							n--;
                                                        //nucs[n] = nucs[--n];
                                                }
                                                //printf("AFTER ALOOP\n");
                                                nucs[a] = singlenuc;
                                                //printf("NUCS: %s", nucs);
                                                size++;
					 }
                                }
                        }
			//printf("NUCS: %s\n", nucs);
                       int low=0, up=NUMNUCS-1, mid;
                       //printf("NUMNUCS: %i", NUMNUCS);
			//printf("Binary search. NUCS: %s\n", nucs);
                        while(1){
                                mid = (up + low)/2;
                                //printf("MID IS: %d\n", mid);
                                //printf("Comparing %s and %s\n", simNucs[mid], nucs);
                                if( strcmp(simNucs[mid], nucs) == 0 ){
                                        index = mid;
                                        break;
                                }
                                else if(up < low){
                                        index = -1;
                                        break;
                                }
                                else if(strcmp(simNucs[mid], nucs) < 0)
                                        low = mid + 1;
                                else
                                        up = mid - 1;
                        }
//printf("Post binary search.\n");
//printf("INDEX: %d\n", index);
//printf("COUNT: %i", count);


			if (index == -1)
				codonids[k] = 15;
			else
				codonids[k] = index;
			if (k%3 == 2)
			{
				if (strcmp(nucs, "A") == 0 || strcmp(nucs, "C") == 0 || strcmp(nucs, "T") == 0 || strcmp(nucs, "G") == 0)
				{	count2++;
//					printf("count2: %d\n", count2);
				}
				if (count2 == 2 || count2 == 3)
				{
					count++;
					numprim++;
//					printf("numprim increment: %d\n", numprim);
				}
				else
				{
					if(count < min_primer_length)
					{
						for (j=0; j<count*3; j++)
							codonids[j+start] = -1;
						numprim -= count;
					}
					start = k+1;
					count = 0;
					codonids[k-1] = -1;
					codonids[k] = -1;
					codonids[k-2] = -1;
				}
				count2 = 0;
			}
			else
			{
				if (strcmp(nucs, "A") == 0 || strcmp(nucs, "C") == 0 || strcmp(nucs, "T") == 0 || strcmp(nucs, "G") == 0)
				{	count2++;
//					printf("count2 in else: %d", count2);
				}
			}

		} // end of k

/*			if(index != -1){
                                codonids[k]=index;
				count++;
				numprim++;
			}

			else{
				if(count < min_primer_length*3){
					for(j=0; j<count; j++)
						codonids[j+start]=-1;
					numprim -= count;
				}
				start=k+1;
				count=0;
			}

*/
//		}	// end of k
		//printf("end of k\n");
		if(count < min_primer_length){
			//printf("COUNT IS: %d\n", count);
			for(j=0; j<count*3; j++)
				codonids[j+start]=-1;
			numprim -= count;
			int mm;
			for (mm=(count*3)+start; mm < seq_length; mm++)
				codonids[mm] = -1; 
		}


		int endseq=0;
		int primercount=0;
//		printf("numprim: %d", numprim);
		if (numprim != 0) {
		printf("=======================================================================================\n");
                printf("PRIMERS FOR GROUP %d:\n", i+groupsb4);
		for(k=0; k < seq_length; k++){
			if(codonids[k] != -1){
				if(endseq==0){
					printf("Primer %i @ %i:", primercount,k);
					primercount++;
					endseq=1;
				}
				printf("%s", nucCodons[codonids[k]]);
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
        } //end of nuc primer creation
}


