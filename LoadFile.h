#ifndef LOADFILE_H
#define LOADFILE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// call LoadFile in cluster.c

void LoadFile(char* filename, char *aligned, FILE *data, int rows, int seq_length);

#endif
