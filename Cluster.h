#ifndef CLUSTER_H
#define CLUSTER_H
//#include "Dimensions.h"

#define GRIDX 65535
#define BLOCKSIZE 1024
#define BLOCKX 1024
#define BLOCKY 1
#define ACIDSIZE 25
#define TABLESIZE 625
//Kernel

__global__ void computeNucSim(bool *similarity_matrix, const char* __restrict sequences, const int* __restrict groups, int N, int M, int P){
//__global__ void computeNucSim(bool *similarity, char *strands, int *groups, int rows, int seq_length){
        // Compute a unique thread number for this thread
   	//unsigned int threadNum = blockIdx.y*GRIDX*BLOCKSIZE + blockIdx.x*BLOCKSIZE + threadIdx.x;
   	unsigned int threadNum = blockIdx.z*GRIDX*GRIDX*BLOCKSIZE + blockIdx.y*GRIDX*BLOCKSIZE + blockIdx.x*BLOCKSIZE + threadIdx.x;

        unsigned int MP = M*P;
   	if(threadNum>=P*MP)
     		return;

        // Map the thread number to (i,j,k) in PxPxM space
    	unsigned int i = threadNum/MP;
        unsigned int offset = threadNum - MP*i;
    	unsigned int j = offset / M;
        if (i > j)
           return;
        unsigned int k = offset % M;
 
        // Find group boundaries for groups i and j
	unsigned int lowerI = (i == 0) ? 0 : groups[i-1]+1;
	unsigned int lowerJ = (j == 0) ? 0 : groups[j-1]+1;
        unsigned int upperI = groups[i];
	unsigned int upperJ = groups[j];

        // Loop over the sequences in each group
        // If we find a single non-similar amino acid at any point
        // set the similarity matrix to zero and return
        unsigned int index = i*M*N+j*M+k;
   	for(unsigned int l=lowerI; l<=upperI; l++){
      		for(unsigned int m=lowerJ; m<=upperJ; m++){
    			char aa1 = sequences[l*M+k];
    			char aa2 = sequences[m*M+k];
			if (aa1 < 'A' || aa2 < 'A' || aa1 != aa2) {
                                        similarity_matrix[index]=0;
                                        return;
                        }
    		}
  	}

        // If we made it all the way through, set the similarity matrix to one
        similarity_matrix[index]=1;
   	/*unsigned int elementNum = blockIdx.y*GRIDX*BLOCKSIZE + blockIdx.x*BLOCKSIZE + threadIdx.x;
   	//unsigned int threadNum = blockIdx.x*GRIDX*BLOCKSIZE + blockIdx.y*BLOCKSIZE + threadIdx.z*BLOCKX*BLOCKY + threadIdx.y*BLOCKX + threadIdx.x;
        //long elementNum = blockIdx.x*65535*BLOCKSIZE + blockIdx.y*BLOCKSIZE + threadIdx.x*64 + threadIdx.y*8 + threadIdx.z;
        //unsigned int elementNum = blockIdx.x*65535*BLOCKSIZE + blockIdx.y*BLOCKSIZE + threadIdx.x*BLOCKX*BLOCKY + threadIdx.y*BLOCKY + threadIdx.z;
        if(elementNum>=(long)rows*rows*seq_length)
                return;
        unsigned int i = elementNum/(rows*seq_length);
        unsigned int j = (elementNum - rows*seq_length*i) / seq_length;
	unsigned int k = (elementNum - rows*seq_length*i) % seq_length;
	
        if (groups[i] == -1 || groups[j] == -1)
		return;
        //long lowerL, lowerM, l, m;
        int lowerL, lowerM, l, m;
        if(i==0)
                lowerL=0;
        else
                lowerL=groups[i-1]+1;
        if(j==0)
                lowerM=0;
        else
                lowerM=groups[j-1]+1;
        for(l=lowerL; l<=groups[i]; l++){
                for(m=lowerM; m<=groups[j]; m++){
                        char aa1 = strands[l*seq_length+k];
                        char aa2 = strands[m*seq_length+k];
                        if (aa1 > aa2){
                                char tmp = aa1;
                                aa1 = aa2;
                                aa2 = tmp;
                        }
                        if (aa1 == 'A'){
                                if (!(aa2 == 'A' )){
                                        similarity[i*rows*seq_length+j*seq_length+k]=0;
                                        return;
                                }
                        }
                        else if(aa1 == '.' || aa1 == '#' || aa1 == '-'){
                                similarity[i*rows*seq_length+j*seq_length+k]=0;
                                return;
			}
			else if(aa1 == 'G'){
		                  if(!(aa2 == 'G')){
			          	similarity[i*rows*seq_length+j*seq_length+k]=0;
				  	return;
				 }
			}
		        else if (aa1 == 'C'){
		        	if(!(aa2 == 'C')){
					similarity[i*rows*seq_length+j*seq_length+k]=0;
					return;
				}
			}
			else if (aa1 == 'T'){
				 if(!(aa2 == 'T')){
			   	        similarity[i*rows*seq_length+j*seq_length+k]=0;
					return;
				}
			}

		} //first for
	} //second for
	similarity[i*rows*seq_length+j*seq_length+k]=1;*/
}

__constant__  unsigned int acidtable[25*25] =
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
       /*A*/ {1,0,0,1,1,0,1,0,0,0,0,1,0,0,0,1,0,0,1,1,0,1,1,0,0, // A
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
       /*B*/  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // B
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
       /*C*/  0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,1, // C
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
       /*D*/  1,0,0,1,1,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1, // D
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
       /*E*/  1,0,0,1,1,0,1,0,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0, // E
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
       /*F*/  0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,1,0,0,1, // F
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
       /*G*/  1,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0, // G
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
       /*H*/  0,0,0,1,0,0,0,1,0,0,0,1,0,1,0,1,1,1,0,0,0,0,0,0,1, // H
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
       /*I*/  0,0,0,0,0,1,0,0,1,0,1,1,1,1,0,0,0,1,1,1,0,1,0,0,0, // I
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
       /*J*/  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // J
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
       /*K*/  0,0,0,0,1,0,0,0,1,0,1,0,1,1,0,0,1,1,0,1,0,0,1,0,0, // K
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
       /*L*/  1,0,0,0,0,1,0,1,1,0,0,1,1,0,0,1,1,1,1,1,0,1,1,0,0, // L
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
       /*M*/  0,0,0,0,0,0,0,0,1,0,1,1,1,0,0,0,0,0,0,1,0,1,1,0,0, // M
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
       /*N*/  0,0,0,1,0,0,0,1,1,0,1,0,0,1,0,0,0,0,1,1,0,0,0,0,1, // N
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
       /*O*/  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // O
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
       /*P*/  1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,1,1,1,1,0,0,1,0,0, // P
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
       /*Q*/  0,0,0,0,1,0,0,1,0,0,1,1,0,0,0,1,1,1,0,0,0,0,1,0,0, // Q
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
       /*R*/  0,0,1,0,0,0,1,1,1,0,1,1,0,0,0,1,1,1,1,1,0,0,1,0,0, // R
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
       /*S*/  1,0,1,0,0,1,1,0,1,0,0,1,0,1,0,1,0,1,1,1,0,0,1,0,1, // S
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
       /*T*/  1,0,0,0,0,0,0,0,1,0,1,1,1,1,0,1,0,1,1,1,0,0,1,0,0, // T
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
       /*U*/  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0, // U
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
       /*V*/  1,0,0,1,1,1,1,0,1,0,0,1,1,0,0,0,0,0,0,0,0,1,1,0,0, // V
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
       /*W*/  1,0,1,0,1,0,1,0,0,0,1,1,1,0,0,1,1,1,1,1,0,1,1,0,0, // W
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
       /*X*/  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // X
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
       /*Y*/  0,0,1,1,0,1,0,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1}; // Y
            //A B C D E F G H I J K L M N O P Q R S T U V W X Y
//#define I(X) (X-65)

        /*__shared__ unsigned int s_acidtable[TABLESIZE];
        if (acidNum < TABLESIZE) {
            s_acidtable[acidNum] = acidtable[acidNum];
        }
        __syncthreads();*/

__global__ void computeSim(bool *similarity_matrix, const char* __restrict sequences, const int* __restrict groups, int N, int M, int P){
        // Compute a unique thread number for this thread
   	//unsigned int threadNum = blockIdx.y*GRIDX*BLOCKSIZE + blockIdx.x*BLOCKSIZE + threadIdx.x;
   	unsigned int threadNum = blockIdx.z*GRIDX*GRIDX*BLOCKSIZE + blockIdx.y*GRIDX*BLOCKSIZE + blockIdx.x*BLOCKSIZE + threadIdx.x;
        unsigned int MP = M*P;
   	if(threadNum>=P*MP)
     		return;

        // Map the thread number to (i,j,k) in PxPxM space
    	unsigned int i = threadNum/MP;
        unsigned int offset = threadNum - MP*i;
    	unsigned int j = offset / M;
        if (i > j)
           return;
        unsigned int k = offset % M;
 
        // Find group boundaries for groups i and j
	unsigned int lowerI = (i == 0) ? 0 : groups[i-1]+1;
	unsigned int lowerJ = (j == 0) ? 0 : groups[j-1]+1;
        unsigned int upperI = groups[i];
	unsigned int upperJ = groups[j];
        
        // Loop over the sequences in each group
        // If we find a single non-similar amino acid at any point
        // set the similarity matrix to zero and return
        unsigned int index = i*M*N+j*M+k;
   	for(unsigned int l=lowerI; l<=upperI; l++){
      		for(unsigned int m=lowerJ; m<=upperJ; m++){
    			char aa1 = sequences[l*M+k];
    			char aa2 = sequences[m*M+k];
			if (aa1 < 'A' || aa2 < 'A' || acidtable[(aa1-65)*ACIDSIZE+(aa2-65)] == 0) {
                                        similarity_matrix[index]=0;
                                        return;
                        }
    		}
  	}

        // If we made it all the way through, set the similarity matrix to one
        similarity_matrix[index]=1;
}


__global__ void computeBlock(const bool* __restrict similarity_matrix, int* block_score, int N, int M, int P, int L){
        // Compute a unique thread number for this thread
        unsigned int threadNum= blockIdx.z*GRIDX*GRIDX*BLOCKSIZE + blockIdx.y*GRIDX*BLOCKSIZE + blockIdx.x*BLOCKSIZE + threadIdx.z*BLOCKX*BLOCKY +threadIdx.y*BLOCKX + threadIdx.x;
        if(threadNum >= P*P)
                return;

        // Map the thread number to (i, j) in PxP space
        unsigned int i = threadNum / P;
        unsigned int j = threadNum % P;
        if (i >= j)
           return;

        // Compute starting indices in similarity matrix for:
        // Groups i and j, Group i with itself, and Group j with itself
        unsigned int IM = i*M;
        unsigned int JM = j*M;
        unsigned int INM = IM*N;
        unsigned int JNM = JM*N;
        unsigned int indexIJ = INM+JM;
        unsigned int indexII = INM+IM;
        unsigned int indexJJ = JNM+JM;
	
        // Set scores to zero, and index for block score of groups i and j
        int score = 0, final_score = 0;
        unsigned int index = i*P+j;

        // Loop over entries of similarity matrix for each group
        for(unsigned int k=0; k<M; k++){
                // If amino acid k is similar when comparing all sequences in i and j,
                // all sequences in i and all sequences in j count it towards the score
                if(similarity_matrix[indexIJ+k] == 1 &&
                   similarity_matrix[indexII+k] == 1 &&
                   similarity_matrix[indexJJ+k] == 1)
                        score++;
                // Otherwise check if our score is large enough, if so count it
                // in the overall score and reset
                else {
			if (score >= L)
                           final_score += score;
                        score = 0;
                }
        }

        // Set block score entry to the overall score, accounting for trailing similar acids
        if(score >= L)
               final_score += score;
        block_score[index] = final_score;
}

__global__ void computeNucBlock(const bool* __restrict similarity, int* block_score, int P, int N, int M, int L){
        // Compute a unique thread number for this thread
        unsigned int threadNum= blockIdx.z*GRIDX*GRIDX*BLOCKSIZE + blockIdx.y*GRIDX*BLOCKSIZE + blockIdx.x*BLOCKSIZE + threadIdx.z*BLOCKX*BLOCKY +threadIdx.y*BLOCKX + threadIdx.x;
        if(threadNum >= P*P)
                return;

        // Map the thread number to (i, j) in PxP space
        unsigned int i = threadNum / P;
        unsigned int j = threadNum % P;
        if (i >= j)
           return;

        // Compute starting indices in similarity matrix for:
        // Groups i and j, Group i with itself, and Group j with itself
        unsigned int IM = i*M;
        unsigned int JM = j*M;
        unsigned int INM = IM*N;
        unsigned int JNM = JM*N;
        unsigned int indexIJ = INM+JM;
        unsigned int indexII = INM+IM;
        unsigned int indexJJ = JNM+JM;
	
        // Set scores to zero, and index for block score of groups i and j
        int score = 0, final_score = 0;
        unsigned int index = i*P+j;
        int k;
        int tempScore;
        /*unsigned int elementNum= blockIdx.x*BLOCKSIZE + threadIdx.z*BLOCKX*BLOCKY +threadIdx.y*BLOCKX + threadIdx.x;
        unsigned int i = elementNum / r;
        unsigned int j = elementNum % r;
        int k;
        if(elementNum >= r*r)
                return;
        if(i >= j)
                *(blockScore+r*i+j) = -1;
        else{
                *(blockScore+r*i+j)=0;
                int score = 0, final_score = 0, tempScore;
                //int primer_length = min_primer_length; //60
                unsigned int temp1 = i*seq_length;
                unsigned int temp2 = temp1*rows;
                unsigned int temp3 = j*seq_length;
                unsigned int temp4 = temp3*rows;
                unsigned int index1 = temp2+temp3;
                unsigned int index2 = temp2+temp1;
                unsigned int index3 = temp4+temp3;
*/
		for(k=0; k<((M/3)*3); k=k+3){
			tempScore=0;
         		if(similarity[indexIJ+k] == 1 &&
			   similarity[indexII+k] == 1 &&
			   similarity[indexJJ+k] == 1) {
                        //if(similarity[i*rows*seq_length+j*seq_length+k] == 1){
				tempScore++;
			}
         		if(similarity[indexIJ+k+1] == 1 &&
			   similarity[indexII+k+1] == 1 &&
			   similarity[indexJJ+k+1] == 1) {
			//if(similarity[i*rows*seq_length+j*seq_length+k+1] == 1){
				tempScore++;
			}
         		if(tempScore < 2 || similarity[indexIJ+k+2] == 1 &&
			   similarity[indexII+k+2] == 1 &&
			   similarity[indexJJ+k+2] == 1) {
			//if(similarity[i*rows*seq_length+j*seq_length+k+2] == 1){
				tempScore++;
			}
                        if((tempScore == 2) || (tempScore ==3)){ // if two or three are the same, then add one to the score (one more similar acid
				score++;
			}
                        else{ // gap in alignment, check if score is > min primer length
                                if(score >= L)
				     final_score += score;
                                 score = 0;  // reset to 0 to keep counting
                        }
                }
                // Make sure the block score was set even if 1 in last position
                if(score >= L)
                        final_score += score;;
	block_score[index] = final_score;
//                        *(blockScore+r*i+j) = final_score;
//        }
}

#endif




