#ifndef COMPAREGROUPK_H
#define COMPAREGROUPK_H

#include <stdlib.h>
#include <stdio.h>
#include "LinkedList.h"
#include "pointerMath.h"
#include "intLinkedList.h"

int compareGroupK(int*, LinkedList* , int*, int* , int* , int* , int*, intLinkedList*);
void fuzzyK(int* blockScore, LinkedList *h_aligned, int* r, int* numCents, int* fuzzymatch, int* cents, int* fuzziness, int* fuzzies, intLinkedList*, int* tolerance);

#endif
