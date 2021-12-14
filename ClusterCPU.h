#ifndef CPUCLUSTER_H
#define CPUCLUSTER_H
//#include "Dimensions.h"

int computeNucSim(char** strands, int* groups, int i, int j, int k);

//Kernel
//int computeSim(char strands[ROWS][SEQ_LENGTH], int groups[ROWS], int i, int j, int k)
//int computeSim(char** strands, int* groups, int i, int j, int k);
void computeSimCPU(char* strands, int* groups, int r, int N, int M, int* similarity_matrix);

void computeNucBlock(char** strands, int* groups, int* blockScore, int r, int SEQ_LENGTH, int MIN_PRIMER_LENGTH);

//void computeBlock(char strands[ROWS][SEQ_LENGTH], int groups[ROWS], int blockScore[ROWS][ROWS])
//void computeBlock(char** strands, int* groups, int* blockScore, int r, int SEQ_LENGTH, int MIN_PRIMER_LENGTH);
void computeBlockCPU(/*char* strands, int* groups,*/ int* similarity_matrix,  int* blockScore, int r, int N, int SEQ_LENGTH, int MIN_PRIMER_LENGTH);

#endif
