#ifndef GETNUCPRIMER_H
#define GETNUCPRIMER_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//#include "Dimensions.h"
#include "LinkedList.h"
#include "NucCodons.h"


void getNucPrimer(LinkedList*,int*,char *align,int, char *simNucs[NUMNUCS], char *nucCodons[NUMNUCS], int seq_length, int min_primer_length);

#endif




