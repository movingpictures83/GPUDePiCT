#include "NucCodons.h"

static int nucPosition = 0;

void insertNucCodon(char* nucs, char* nucCodon, char* simNucs[NUMNUCS], char* nucCodons[NUMNUCS])
{
   simNucs[nucPosition] = nucs;
   nucCodons[nucPosition] = nucCodon;
   nucPosition++;
}

void insertNucCodons(char* simNucs[NUMNUCS], char* nucCodons[NUMNUCS])
{
// ambiguity codes from: http://droog.gs.washington.edu/parc/images/iupac.html
   insertNucCodon("A", "T", simNucs, nucCodons);
   insertNucCodon("AC", "K", simNucs, nucCodons);
   insertNucCodon("ACG", "B", simNucs, nucCodons);
   insertNucCodon("ACGT", "N", simNucs, nucCodons);
   insertNucCodon("ACT", "D", simNucs, nucCodons);
   insertNucCodon("AG", "Y", simNucs, nucCodons);
   insertNucCodon("AGT", "H", simNucs, nucCodons);
   insertNucCodon("AT", "W", simNucs, nucCodons);
   insertNucCodon("C", "G", simNucs, nucCodons);
   insertNucCodon("CG", "S", simNucs, nucCodons);
   insertNucCodon("CGT", "V", simNucs, nucCodons);
   insertNucCodon("CT", "R", simNucs, nucCodons);
   insertNucCodon("G", "C", simNucs, nucCodons);
   insertNucCodon("GT", "M", simNucs, nucCodons);
   insertNucCodon("T", "A", simNucs, nucCodons);
   insertNucCodon("dum", "N", simNucs, nucCodons);

}
