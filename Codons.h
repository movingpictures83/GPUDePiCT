#ifndef CODONS_H
#define CODONS_H

#define NUMACIDS 309

//static char* simacids[NUMACIDS];
//static char* codons[NUMACIDS];
static int position = 0;
void insertCodon(char* acids, char* codon, char* simacids[NUMACIDS], char* codons[NUMACIDS]);

void insertCodons(char* simacids[NUMACIDS], char* codons[NUMACIDS]);

#endif
