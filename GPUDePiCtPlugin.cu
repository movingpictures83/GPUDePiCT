#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "LinkedList.h"
#include "intLinkedList.h"
#include <sys/time.h>
#include "ClusterCPU.h"
#include "Codons.h"
#include "LoadFile.h"
#include "NucCodons.h"
#include "flatten.h"
#include "compareGroup.h"
#include "compareGroupK.h"
#include "getprimer.h"
#include "getNucPrimer.h"
#include "pointerMath.h"
#include "Cluster.h"
/*  Degenerate Primer Design via Clustering,		  *
 *  on a GPU utilizing Nvidia's Cuda Platform		  *
 *  		  					  *
 *  Written by Philippe Novikov & James Parda		  *
 *  Summer Research 2012 with Dr. Trevor Cickovski        *
 *                                                        */
#include "GPUDePiCtPlugin.h"
#include <fstream>

void GPUDePiCtPlugin::input(std::string file) {
 inputfile = file;
 std::ifstream ifile(inputfile.c_str(), std::ios::in);
 while (!ifile.eof()) {
   std::string key, value;
   ifile >> key;
   ifile >> value;

   parameters[key] = value;
 }
 std::string composite = std::string(PluginManager::prefix())+"/"+parameters["inputfile"];
		filename = (char*) composite.c_str();
		ROWS = atoi(parameters["rows"].c_str());
		SEQ_LENGTH = atoi(parameters["seqlength"].c_str());
		boolVal = atoi(parameters["AAorNuc"].c_str());
		algorithm = atoi(parameters["kmeans"].c_str());
		fuzziness = atoi(parameters["fuzzy"].c_str());
		tolerance = atoi(parameters["tolerance"].c_str());
		if (parameters.count("minlength") != 0)
			MIN_PRIMER_LENGTH = atoi(parameters["minlength"].c_str());

	data = fopen(filename, "r");
	aligned = (char*) malloc(ROWS*SEQ_LENGTH*sizeof(char));

	// filetype checker and reader
	LoadFile(filename, aligned, data, ROWS, SEQ_LENGTH);
}

void GPUDePiCtPlugin::run() {
	//load correct set of codons
  	if(boolVal==0)
		insertCodons(simacids, codons);
	else
		insertNucCodons(simNucs, nucCodons);

	r = ROWS;
	int *done = (int*)malloc(ROWS*sizeof(int));


/* Cluster sequences */
  	//Start timer
  	struct timeval start, finish;
  	struct timezone tz;
  	gettimeofday(&start, &tz);

      	// Populate the list with our sequences
      	h_aligned = (LinkedList*)malloc(sizeof(LinkedList));
      	Node *root = (Node*)malloc(sizeof(Node));
      	CharNode *first_sequence = (CharNode*)malloc(sizeof(CharNode));
      	first_sequence->data = &aligned[0*SEQ_LENGTH+0];
	first_sequence->next = NULL;
      	root->data = first_sequence;
        root->next = NULL;
	intLinkedList *tracer = (intLinkedList*)malloc(sizeof(intLinkedList));
	tracer->root = NULL;
	for(x=0; x<ROWS; x++){
		intNodeI *xintnode=(intNodeI*)malloc(sizeof(intNodeI));
		xintnode->data=x;
		xintnode->next=NULL;
		NodeI *xnodeI =(NodeI*)malloc(sizeof(NodeI));
		xnodeI->data=xintnode;
		xnodeI->next=NULL;
		insertendI(tracer, xnodeI);
	}
	//List_printI(tracer);
	intLinkedList *final_tracer =(intLinkedList*)malloc(sizeof(intLinkedList));
        final_tracer->root = NULL;
	h_aligned->root = root;
        h_aligned->num_elems = 1;
        for(x=1; x<ROWS; x++) {
		int q;
		/*printf("INSERTING STRING: ");
		for (q = 0; q < SEQ_LENGTH; q++)
		   printf("%c", aligned[x*SEQ_LENGTH+q]);
		printf("\n");*/
                insert(h_aligned, x+1, &aligned[x*SEQ_LENGTH+0]);
        }
	groups = (int*) malloc(ROWS*sizeof(int));

	//initialize second list for when groups with a max similarity of 0 are removed
	final_list = (LinkedList*)malloc(sizeof(LinkedList));
	final_list->root = NULL;
	final_groups = (int*) malloc(ROWS*sizeof(int));
	for (i = 0; i < ROWS; i++)
		final_groups[i] = 0;
	final_list->num_elems=0;

	#ifndef CPU
		//declaration of gpu data structures
    		cudaMalloc( (void**) &gpu_sim, (unsigned int)ROWS*ROWS*SEQ_LENGTH*sizeof(bool));
    		cudaMalloc( (void**) &gpu_group, ROWS*sizeof(int));
    		//size_t size = SEQ_LENGTH*sizeof(char); // size of an alignment
		cudaMalloc( (void**) &gpu_aligned, ROWS*SEQ_LENGTH*sizeof(char));

		//Define numBlocks and block dimensions for gpu operations
      		unsigned int numBlocks = ((unsigned int)ROWS*ROWS*SEQ_LENGTH - 1)/BLOCKSIZE+1;
      		unsigned int numBlocks2 =((unsigned int)ROWS*ROWS-1)/BLOCKSIZE+1;
		unsigned int numGrids=1;
		unsigned int numGridsY=1;
		unsigned int numGrids2=1;
		unsigned int numGrids2Y=1;
      		if (numBlocks > 65535){
         		numGrids = numBlocks / 65535;
         		if (numGrids % numBlocks != 0)
				numGrids++;
         		numBlocks = 65535;
			if (numGrids > 65535) {
				numGridsY = numGrids / 65535;
				if (numGridsY % numGrids != 0)
					numGridsY++;
				numGrids = 65535;
			}
      		}
		if (numBlocks2 > 65535) {
			numGrids2 = numBlocks2 / 65535;
			if (numGrids2 % numBlocks2 != 0)
				numGrids2++;
			numBlocks2 = 65535;
			if (numGrids2 > 65535) {
				numGrids2Y = numGrids2 / 65535;
				if (numGrids2Y % numGrids2 != 0)
					numGrids2Y++;
				numGrids2 = 65535;
			}
		}
		int BLOCKXY = BLOCKX*BLOCKY;
		dim3 dimBlock(BLOCKX, BLOCKY, BLOCKSIZE/BLOCKXY);
		//dim3 dimBlock(BLOCKSIZE/BLOCKXY, BLOCKX, BLOCKY);
		dim3 dimBlock2(BLOCKX, BLOCKY, BLOCKSIZE/BLOCKXY);
		//dim3 dimBlock2(BLOCKSIZE/BLOCKXY, BLOCKX, BLOCKY);
		dim3 dimGrid(65535, numGrids, numGridsY);
		dim3 dimGrid2(65535, numGrids2, numGrids2Y);
      		printf("%d grids of %d blocks have been allocated for this process.\n", numGrids, numBlocks);
		printf("%d blocks have been allocated for this process.\n",numBlocks2);
	#else
	   int* sim_mat = malloc((int)ROWS*ROWS*SEQ_LENGTH*sizeof(int));
	#endif
			/*printf("**********************************************************************************\n");
			printf("* OUR LIST IN MAIN                                                                \n");
			List_print(h_aligned, 1126);
			printf("**********************************************************************************\n");*/
	char* cpu_aligned = (char*)malloc(ROWS*SEQ_LENGTH*sizeof(char));
	//char cpu_aligned[ROWS][SEQ_LENGTH];
    	while(1){
		//flatten h_aligned and copy it to the gpu
		/*int a, c;;
		for (a = 0; a < ROWS; a++) {
			for (c = 0; c < SEQ_LENGTH; c++) {
				printf("%c ", aligned[a*SEQ_LENGTH+c]);
			}
			printf("\n");
		}*/
      		flatten(h_aligned->root, cpu_aligned, groups, ROWS, SEQ_LENGTH);
		#ifndef CPU
 	//			bool* cpu_sim = (bool*)malloc(ROWS*ROWS*SEQ_LENGTH*sizeof(bool));

			cudaMemcpy(gpu_aligned, cpu_aligned, ROWS*SEQ_LENGTH*sizeof(char), cudaMemcpyHostToDevice);
      			cudaMemcpy(gpu_group, groups, ROWS*sizeof(int), cudaMemcpyHostToDevice);
			//compute sim for either aa or nuc
			if(boolVal ==0)
				computeSim<<< dimGrid, dimBlock >>>( (bool (*))gpu_sim, (char (*))gpu_aligned, (int (*))gpu_group, ROWS, SEQ_LENGTH, r );
			else
				computeNucSim<<< dimGrid, dimBlock >>>( (bool (*))gpu_sim, (char (*))gpu_aligned, (int (*))gpu_group, ROWS, SEQ_LENGTH, r );
		//initialize blockScore(cpu) and bScore(gpu)

/*		cudaMemcpy(cpu_sim, gpu_sim, ROWS*ROWS*SEQ_LENGTH*sizeof(bool), cudaMemcpyDeviceToHost);
		int j, k;
				for (k=0; k<SEQ_LENGTH; k++){
					printf("%d", cpu_sim[0*ROWS*SEQ_LENGTH+1*SEQ_LENGTH+k]);
				}
				printf("\n");
				List_print(h_aligned, SEQ_LENGTH);
*/		
		#else
			//printf("COMPUTING SIM....\n");
			computeSimCPU(cpu_aligned, groups, r, ROWS, SEQ_LENGTH, sim_mat);
			/*int a, b, c;
			for (a = 0; a < ROWS; a++)
			for (b = 0; b < ROWS; b++)
			for (c = 0; c < SEQ_LENGTH; c++)
			   if (sim_mat[a*ROWS*SEQ_LENGTH+b*SEQ_LENGTH+c] != 0)
			      printf("NON-ZERO ENTRY FOUND\n");*/
			//printf("DONE....\n");
		#endif

		int* blockScore=(int*)malloc(r*r*sizeof(int));
		#ifndef CPU
			cudaMalloc((void**)&bScore, r*r*sizeof(int));
			//printf("BEFORE BLOCK\n");
			//List_print(h_aligned, SEQ_LENGTH);
			//compute block score for nuc or aa, copy back to cpu and free gpu version
			//printf("ALLOCATED %d\n", r*r);
			if(boolVal ==0)
				computeBlock<<<dimGrid2, dimBlock2>>>((bool (*))gpu_sim, (int  (*))bScore, ROWS, SEQ_LENGTH, r, MIN_PRIMER_LENGTH);
      			else{
				//printf("COMPUTE NUC BLOCK %d\n", numBlocks2);
				computeNucBlock<<<dimGrid2, dimBlock2>>>((bool (*))gpu_sim, (int  (*))bScore, r, ROWS, SEQ_LENGTH, MIN_PRIMER_LENGTH);
      			}
			/*gettimeofday(&finish, &tz);
        		double block1Elapsed = ( finish.tv_sec - start.tv_sec ) * 1000.0 + ( finish.tv_usec - start.tv_usec ) / 1000.0;
        		printf("Time to before copy: %lf ms\n", block1Elapsed);*/

			cudaMemcpy(blockScore, bScore, r*r*sizeof(int), cudaMemcpyDeviceToHost);
			cudaFree(bScore);
		#else
			//printf("R:%d\n", r);
			if(boolVal ==0)
				computeBlockCPU(sim_mat,/*cpu_aligned, groups,*/ blockScore, r, ROWS, SEQ_LENGTH, MIN_PRIMER_LENGTH);
			else
				computeNucBlockCPU(cpu_aligned, groups, blockScore, r, SEQ_LENGTH, MIN_PRIMER_LENGTH);
		#endif

		/*gettimeofday(&finish, &tz);
        	double block2Elapsed = ( finish.tv_sec - start.tv_sec ) * 1000.0 + ( finish.tv_usec - start.tv_usec ) / 1000.0;
        	printf("Time to after copy: %lf ms\n", block2Elapsed);*/
		/*int v, b;
		printf("    ");
		for(v=0; v<r; v++)
			printf("G%2d|", v);
		printf("\n");
		for(v=0; v<r; v++){
			printf("G%2d|", v);
			for(b=0; b<r; b++)
				printf("%3d|", *(blockScore+v*r+b));
			printf("\n");
		}*/


		//initialize arrays for centroids
		//int done[r];
		int moreCents[r];
		//int i;
		for(i=0; i<r; i++){
			moreCents[i]=-1;
			//done[i]=0;
		}
		numCents=0;
		int cents[r];
		for(i=0; i<r; i++){
			cents[i*2+0]=-1;
			cents[i*2+1]=-1;
		}

		//int *done = (int*)malloc(ROWS*sizeof(int));
		for (i = 0; i < ROWS; i++)
			done[i] = 0;

		//initializes arrays for fuzzy K
		int fuzzymatch[r*fuzziness];
		for(i=0; i<r*fuzziness; i++)
			fuzzymatch[i]=-1;
		fuzzies=0;
		//finds centroids and matches and merges them, returns 1 if still more groups to merge, -1 if not
		//int max_j;
		if (algorithm == 0)
			max_j=compareGroup(blockScore, h_aligned, done, &r, tracer);
		else{
			max_j=compareGroupK(blockScore, h_aligned, done, &r, moreCents, cents, &numCents, tracer);
			if(fuzziness!=0&&max_j!=-1)
				fuzzyK(blockScore, h_aligned, &r, &numCents, fuzzymatch, cents, &fuzziness, &fuzzies, tracer, &tolerance);
    		}
		//printf("MAX J: %d", max_j);
		if(max_j==-1){
			//printf("max_j -1");
			break;
		}
		//printf("after compare group\n");
		//#ifdef CPU
			free(blockScore);
		//#endif
		//continue
		if (algorithm == 1) {
			printf("Number of centroids:%d\n", numCents);
				/*for(i=0; i<numCents; i++)
				printf("Centroid %d: %d Match: %d\n", i, cents[i*2+0], cents[i*2+1]);*/
			printf("Number of fuzzy matches:%d\n", fuzzies);
		}
		//remove groups similar to no other groups
    		//int count;
		w=0;
		//int stop;
		if (algorithm==0) stop = r+1;
		else stop = r;
		for(count=0; count<stop; count++){
			//printf("count: %d\n", count);
			if (algorithm == 1) {
				for(i=0; i<numCents; i++){
					//accounts for centroid matches already removed from list
					if(count==cents[i*2+1]){
						//printf("centmatch: %d\n", cents[i*2+1]);
						w++;
					}
					//int f;
					for(f=0; f<fuzziness; f++){
						if(count==fuzzymatch[i*fuzziness+f]){
							//printf("fuzzymatch: %d\n", fuzzymatch[i*fuzziness+f]);
							w++;
						}
						else if(fuzzymatch[i*fuzziness+f]==-1)
							break;
					}
				}
			}
			else {
				 if (count == max_j) {
					w++;
				}
			}
			//remove group from h_aligned, and add to final_list and update final_groups
			//printf("DONE[%d] = %d\n", count, done[count]);
			if(done[count]==0){
				Node* bad;
				NodeI* badT;
				printf("Removed:(%d)%d\n", count, count-w);
				bad=rmgetnode(h_aligned, count-w);
				badT=rmgetnodeI(tracer, count-w);
				//int groupsize;
				if (count == 0) groupsize = groups[count]+1;
				else groupsize = groups[count]-groups[count-1];
				if(final_list->num_elems > 0)
					final_groups[final_list->num_elems]=final_groups[(final_list->num_elems)-1]+groupsize;
				else
					final_groups[final_list->num_elems]=groupsize-1;
				(final_list->num_elems)++;
				insertend(final_list, bad);
				insertendI(final_tracer, badT);
				free(bad);
				free(badT);
				w++;
				h_aligned->num_elems--;
				tracer->num_elems--;
			}
		}
		//keep r current
		r=h_aligned->num_elems;
		//printf("AFTER BLOCK\n");
		//List_print(h_aligned, SEQ_LENGTH);

		/*printf("Tracer:\n");
		List_printI(tracer);
		printf("Final Tracer:\n");
		List_printI(final_tracer);*/
  	}
     	fclose(data);

	//stop timer
	gettimeofday(&finish, &tz);
     	double clusterElapsed = ( finish.tv_sec - start.tv_sec ) * 1000.0 + ( finish.tv_usec - start.tv_usec ) / 1000.0;
     	//printf("Time to cluster: %lf ms\n", clusterElapsed);
	//printf("Tracer\n");
	//List_printI(tracer);
	//printf("Final tracer\n");
	//List_printI(final_tracer);
/* Design Primers */
 	printf("Designing primers...\n");
}

void GPUDePiCtPlugin::output(std::string file) {
	// flatten final_list and h_aligned
	char* final_aligned = (char*)malloc(ROWS*SEQ_LENGTH*sizeof(char));;
	flatten(final_list->root, final_aligned, NULL, ROWS, SEQ_LENGTH);
      	flatten(h_aligned->root,  aligned, groups, ROWS, SEQ_LENGTH);

	//generate primers for both lists
	if(boolVal==0){
		getprimer(h_aligned, groups, aligned, 0, simacids, codons, SEQ_LENGTH, MIN_PRIMER_LENGTH);
		getprimer(final_list, final_groups, final_aligned, h_aligned->num_elems, simacids, codons, SEQ_LENGTH, MIN_PRIMER_LENGTH);
     	}
	else{
		getNucPrimer(h_aligned, groups, aligned, 0, simNucs, nucCodons, SEQ_LENGTH, MIN_PRIMER_LENGTH);
        	getNucPrimer(final_list, final_groups, final_aligned, h_aligned->num_elems, simNucs, nucCodons, SEQ_LENGTH, MIN_PRIMER_LENGTH);
	}

	//free data
//	free(aligned);
//	free(done);
//	free_list(h_aligned);
//	free_Ilist(tracer);
//	free_Ilist(final_tracer);
//	free(groups);
//	free_list(final_list);
//	free(final_groups);
//	free(cpu_aligned);
//	free(final_aligned);

	#ifndef CPU
		cudaFree(gpu_aligned);
		cudaFree(gpu_sim);
     		cudaFree(gpu_group);
	#endif
}

PluginProxy<GPUDePiCtPlugin> GPUDePiCtPluginProxy = PluginProxy<GPUDePiCtPlugin>("GPUDePiCt", PluginManager::getInstance());


