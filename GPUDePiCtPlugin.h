#include "Plugin.h"
#include "PluginProxy.h"
#include <string>
#include <map>

class GPUDePiCtPlugin : public Plugin {

	public:
		void input(std::string file);
		void run();
		void output(std::string file);
	private:

        int r, x, i, boolVal, algorithm, fuzziness, tolerance;
        int numCents, fuzzies, max_j, count, w, stop, f, groupsize;

        char* simNucs[NUMNUCS];
        char* nucCodons[NUMNUCS];

        char* simacids[NUMACIDS];
        char* codons[NUMACIDS];

        char* filename;
        int ROWS;
        int SEQ_LENGTH;
        int MIN_PRIMER_LENGTH = 3;
	 std::map<std::string, std::string> parameters;
        std::string inputfile;
	std::string outputfile;
	char* aligned;
	FILE* data;
	LinkedList* final_list;
	LinkedList* h_aligned;
	int* groups;
	int* final_groups;
		bool* gpu_sim;
		int* bScore;
		unsigned int *gpu_group;
		char* gpu_aligned;
};
