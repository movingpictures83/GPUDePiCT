#include "LoadFile.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// call LoadFile in cluster.c

void LoadFile(char* filename, char *aligned, FILE *data, int rows, int seq_length)
{

//printf("got into LoadFile.h\n");
	int length;

	// get length of file name
	length = strlen(filename);
	// check file type

//	printf("length: %i\n", length);
//	printf("length-3: %c\n", filename[length-3]);
//	printf("length-2: %c\n", filename[length-2]);
//	printf("length-1: %c\n", filename[length-1]);


	// txt reader
	if (filename[length-3] == 't' && filename[length-2] == 'x' && filename[length-1] == 't')
	{
		int r, x;
		// load example sequences
		for (x=0; x<rows; x++)
		{
			char newline;
			for (r=0; r<seq_length; r++)
			{
				fscanf(data,"%c", &aligned[x*seq_length+r]);
				aligned[x*seq_length+r] = toupper(aligned[x*seq_length+r]);
			}
			fscanf(data, "%c", &newline);
		}
	}

	// msf reader
	else if (filename[length-3] == 'm' && filename[length-2] == 's' && filename[length-1] == 'f')
	{
		size_t size = 500;
		char* line = (char*)malloc(500*sizeof(char));
		char* name = (char*)malloc(500*sizeof(char));
		while (getline(&line, &size, data) != EOF)
		{
			int counter = 0;
			int i, j, k, b, remaining, x;
			char amino, newline;

//			printf("LINE: %s\n", line);
//			printf("%d %c %c\n", (int)strlen(line), line[0], line[1]);

			// search for '//'
			if (strlen(line) >= 2 && line[0]=='/' && line[1]=='/')
			{
//				printf("DOING IF\n");

				x = seq_length/50;
				// loop through name/strand lines
				for (i=1; i<=x; i++)
				{
					for (j=0; j<rows; j++)
					{
						// get name
						fscanf(data, "%s", name);

						// checks if name is sequence length description
						int name_length = strlen(name);
						int num_count = 0;
						for (b=0; b<name_length; b++)
						{
							int compared = (int)name[b];
							if (isdigit(compared))
							{
								num_count++;
							}
						}
						if (num_count == name_length)
						{
							fscanf(data, "%s", name);
							fscanf(data, "%c", &newline);
							fscanf(data, "%s", name);
						}

						// get amino acids and put in array
						for (k=0; k<50; k++)
						{
							fscanf(data, "%c", &amino);
							while (amino == ' ' || amino == '\t')
							{
								fscanf(data, "%c", &amino);
							}
							aligned[j*seq_length+counter+k] = toupper(amino);
						}
						fscanf(data, "%c", &newline);
					}
					counter += 50;
				}
				// read remaining acids
				remaining = seq_length % 50;
				for (j=0; j<rows; j++)
				{
					fscanf(data, "%s", name);
                                        // checks if name is sequence length description
                                        int name_length = strlen(name);
                                        int num_count = 0;
                                        for (b=0; b<name_length; b++)
                                        {
						int compared = (int)name[b];
						if (isdigit(compared))
                                                {
                                                	num_count++;
                                                }
                                        }
                                        if (num_count == name_length)
                                        {
                                                fscanf(data, "%s", name);
                                                fscanf(data, "%c", &newline);
                                                fscanf(data, "%s", name);
                                        }

					for (k=0; k<remaining; k++)
					{
						fscanf(data, "%c", &amino);
						while (amino == ' ' || amino == '\t')
						{
							fscanf(data, "%c", &amino);
						}
						aligned[j*seq_length+counter+k] = toupper(amino);
					}
					fscanf(data, "%c", &newline);
				}
			}
		}
	}

	// fasta reader
	else if (filename[length-5] == 'f' && filename[length-4] == 'a' && filename[length-3] == 's' && filename[length-2] == 't' && filename[length-1] == 'a')
	{
		int a, b; // x, c, remaining, counter, d;
		char amino, newline; // indent;

		char* name = (char*)malloc(500*sizeof(char));
		//char* line = (char*)malloc(500*sizeof(char));
		size_t size = 500;

		// load rows
		for (a=0; a<rows; a++)
		{

			// get name line
			getline(&name, &size, data);

			// gets lines for each sequence
//			x = seq_length/60;

//			counter = 0;
			// loop through lines of sequence
			for (b=0; b<seq_length; b++)
			{
				// loop through aminos in each line
				//for (c=0; c<60; c++)
				//{
					fscanf(data, "%c", &amino);
					while (amino == ' ' || amino == '\t' || amino == '\n')
					{
						fscanf(data, "%c", &amino);
	//					c++;
					}
					aligned[a*seq_length+b] = toupper(amino);
				}
				fscanf(data, "%c", &newline);
	//			counter+=60;
//			}
/*
			// loop through remaining line
			remaining = seq_length % 60;
			for (d=0; d<remaining; d++)
			{
				fscanf(data, "%c", &amino);
                                        while (amino == ' ' || amino == '\t')
					{
                                                fscanf(data, "%c", &amino);
						c++;
					}
				aligned[a*seq_length+counter+d] = toupper(amino);
			}
			fscanf(data, "%c", &newline);
*/		}

	}

        // clustal-w reader
        else if (filename[length-3] == 'a' && filename[length-2] == 'l' && filename[length-1] == 'n')
	{
		size_t size = 500;
		char* line = (char*)malloc(500*sizeof(char));
		char* name = (char*)malloc(500*sizeof(char));

		int x, j, k, l, counter, remaining;
		char amino;

		// get title line
		getline(&line, &size, data);
		// get sections of sequences
		int numperline = 60; // Default, but we'll know for sure after first read
		x = seq_length/numperline;
		counter = 0;
		// loop through sections
		for (j=0; j<x; j++)
		{
			// loop through each row
			for (k=0; k<rows; k++)
			{
				// get name
				fscanf(data, "%s", name);
				// get spaces up to acid
				fscanf(data, "%c", &amino);
				// ignore white space before acids
				while (amino == ' ' || amino == '\t')
				{
					fscanf(data, "%c", &amino);
				}
				// loop through acids
				numperline = 0; // Don't count newline, may not be 60 (default)
				while (amino != '\n') 
				//for (l=0; l<60; l++)
				{
					aligned[k*seq_length+counter+l] = toupper(amino);
					fscanf(data, "%c", &amino);
					numperline++;
					// note: last amino will be newline but not saved
				}
				// First Time Though
				if (j == 0 && k == 0) {
				   x = seq_length / numperline; // Will affect j loop if numperline is not 60
			           //printf("X IS NOW: %d NUMERPERLINE: %d \n", x, numperline);
				}
			}
			counter += numperline; //60;
		}
                // remaining acids
                remaining = seq_length % numperline;//60;
		// loop through each row
         	for (k=0; k<rows; k++)
         	{
			// get name
        	        fscanf(data, "%s", name);
              		// get spaces up to acid
              		fscanf(data, "%c", &amino);
              		// ignore white space before acids
              		while (amino == ' ' || amino == '\t')
              		{
				fscanf(data, "%c", &amino);
			}
              		// loop through acids
              		for (l=0; l<remaining; l++)
              		{
				aligned[k*seq_length+counter+l] = toupper(amino);
				fscanf(data, "%c", &amino);
				// note: last amino will be newline but not saved
			}
		}
	}

        // stockholm reader
        else if (filename[length-3] == 's' && filename[length-2] == 't' && filename[length-1] == 'o')
        {
		size_t size = 500;
		char* line = (char*)malloc(500*sizeof(char));
		char* name = (char*)malloc(500*sizeof(char));

		int x, j, k, l, counter, remaining;
		char amino;

		// get title line
		getline(&line, &size, data);

		// section count
		x = seq_length/100;
                counter = 0;
                // loop through sections
                for (j=0; j<x; j++)
                {
                        // loop through each row
                        for (k=0; k<rows; k++)
                        {
                                // get name
                                fscanf(data, "%s", name);
				// check if name was commenting line
				while (name[0] == '#')
				{
					getline(&line, &size, data);
					fscanf(data, "%s", name);
				}
                                // get spaces up to acid
                                fscanf(data, "%c", &amino);
                                // ignore white space before acids
                                while (amino == ' ' || amino == '\t')
                                {
                                        fscanf(data, "%c", &amino);
                                }
                                // loop through acids
                                for (l=0; l<100; l++)
                                {
                                        aligned[k*seq_length+counter+l] = toupper(amino);
                                        fscanf(data, "%c", &amino);
                                        // note: last amino will be newline but not saved
                                }
                        }
                        counter += 100;
                }
                // remaining acids
                remaining = seq_length % 100;
                // loop through each row
                for (k=0; k<rows; k++)
                {
                        // get name
                        fscanf(data, "%s", name);
                        // check if name was commenting line
                        while (name[0] == '#')
                        {
                                getline(&line, &size, data);
                                fscanf(data, "%s", name);
                        }
                        // get spaces up to acid
                        fscanf(data, "%c", &amino);
                        // ignore white space before acids
                        while (amino == ' ' || amino == '\t')
                        {
                                fscanf(data, "%c", &amino);
                        }
                        // loop through acids
                        for (l=0; l<remaining; l++)
                        {
                                aligned[k*seq_length+counter+l] = toupper(amino);
                                fscanf(data, "%c", &amino);
                                // note: last amino will be newline but not saved
                        }
                }
        }

	// incorrect file type
	else
	{
	printf("Incorrect file type: must be .txt, .msf, .fasta, .sto\n");
	exit(1);
	}
}
