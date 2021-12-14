#ifndef GETPRIMER_H
#define GETPRIMER_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//#include "Dimensions.h"
#include "LinkedList.h"
#include "Codons.h"

void getprimer(LinkedList*,int*,char *align,int, char *simacids[NUMACIDS], char *codons[NUMACIDS], int seq_length, int min_primer_length);

#endif

