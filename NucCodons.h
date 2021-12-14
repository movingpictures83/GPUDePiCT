#ifndef NUCCODONS_H
#define NUCCODONS_H

#define NUMNUCS 16

//static char* simNucs[NUMNUCS];
//static char* nucCodons[NUMNUCS];
//int nucPosition = 0;

void insertNucCodon(char* nucs, char* nucCodon, char* simNucs[NUMNUCS], char* nucCodons[NUMNUCS]);

void insertNucCodons(char* simNucs[NUMNUCS], char* nucCodons[NUMNUCS]);

#endif
