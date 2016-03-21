// make a color map of sensitivity as a function scintillation time-scale and bandwidth
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include "detection.h"
#include <mpi.h>

int main (int argc, char* argv[])
{
	int id;  //  process rank
	int p;   //  number of processes
	//int num;
	int nMax;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &id);
	MPI_Comm_size (MPI_COMM_WORLD, &p);

	//FILE *fin;
	controlStruct control;
	acfStruct acfStructure;
	noiseStruct noiseStructure;

	char fname[1024];   // read in parameter file
	char Tname[1024];   // read in scintillation timescale
	char Fname[1024];   // read in scintillation bandwidth

	int nt, nf;
	double *tdiss, *fdiss;

	//char oname[1024];   // output file name
	//char pname[1024];   // file to plot
	char dname[1024];   // graphics device

	int i, j;
	int n = 1;  // number of dynamic spectrum to simulate

	double flux0;
	double flux1;
	double tdiff;
	double fdiff;

	int plotMode = 0;  // plot existing dynamic spectrum, off by default
	control.noplot = 1;  // don't show dynamic spectrum while simulationg, on by default

	// read options
	for (i=0;i<argc;i++)
	{
		if (strcmp(argv[i],"-f") == 0)
		{
			strcpy(fname,argv[++i]);
			strcpy(Tname,argv[++i]);
			strcpy(Fname,argv[++i]);
			printf ("Parameters are in %s\n", fname);
		}
		//else if (strcmp(argv[i],"-o")==0)
		//{
		//	strcpy(oname,argv[++i]);
		//	//printf ("Dynamic spectrum is output into %s\n", oname);
		//}
		else if (strcmp(argv[i],"-n")==0)
		{
			n = atoi (argv[++i]);
			printf ("Number of dynamic spectrum to simulate %d\n", n);
		}
		//else if (strcmp(argv[i],"-p")==0) // just plot 
		//{
		//	strcpy(pname,argv[++i]);
		//	plotMode = 1;
		//	printf ("Plotting exiting dynamic spectrum.\n");
		//}
		//else if (strcmp(argv[i],"-noplot")==0)
		//{
		//	control.noplot = 1; // Don't show dynamic spetrum while simulationg
		//}
		//else if (strcmp(argv[i],"-dev")==0)
		//{
		//	strcpy(dname,argv[++i]);
		//	printf ("Graphic device: %s\n", dname);
		//}
	}

	nt = readDissNum (Tname);
	nf = readDissNum (Fname);
	tdiss = (double *)malloc(sizeof(double)*nt);
	fdiss = (double *)malloc(sizeof(double)*nf);
	readDiss (Tname, Fname, tdiss, fdiss);

	sprintf(dname, "%s", "1/xs");
	/*
	if ((fin=fopen(oname, "w"))==NULL)
	{
		printf ("Can't open file...\n");
		exit(1);
	}
	*/

	if (plotMode==0)
	{
		// Simulate dynamic spectrum
		// read parameters
		initialiseControl(&control);
		readParams (fname,dname,n,&control);
		//readParams (fname,oname,dname,n,&control);
		printf ("Finished reading parameters.\n");

		//preAllocateMemory (&acfStructure, &control);
		//allocateMemory (&acfStructure);

		calNoise (&noiseStructure, &control);

		//num = (int)((control.scint_ts1-control.scint_ts0)/control.scint_ts_step);
		//for (tdiff=control.scint_ts0; tdiff<control.scint_ts1; tdiff+=control.scint_ts_step)
		//for (i=id; i<=num; i+=p)
		for (i=id; i<nt; i+=p)
		{
			tdiff = pow(10.0, tdiss[i]);   // log scale
			//tdiff = tdiss[i];
			//tdiff = control.scint_ts0+i*control.scint_ts_step;

			for (j=0; j<nf; j++)
			//for (fdiff=control.scint_freqbw0; fdiff<control.scint_freqbw1; fdiff+=control.scint_f_step)
			{
				fdiff = pow(10.0, fdiss[j]);   // log scale
				//fdiff = fdiss[j];
				//printf ("tdiff fdiff: %lf %lf\n", tdiff, fdiff);
				flux0 = control.cFlux0;
				flux1 = control.cFlux1;
				
				control.scint_ts = tdiff;
				control.scint_freqbw = fdiff;

				acfStructure.probability = 1.0;

				nMax = 0;
				while (fabs(acfStructure.probability-0.8)>=0.01 && nMax <= 100)
				{
					control.cFlux = flux0+(flux1-flux0)/2.0;
					//printf ("%lf %lf %.8lf %.3f\n", tdiff, fdiff, control.cFlux, acfStructure.probability);
		
					// simulate dynamic spectra
					calculateScintScale (&acfStructure, &control);

					//if (control.noplot==0 && n == 1)
					//{
					//	// plot while simulating
					//	heatMap (&acfStructure, dname);
					//}

					qualifyVar (&acfStructure, &noiseStructure, &control);

					if (acfStructure.probability>0.8)
					{
						flux1 = control.cFlux;
					}
					else 
					{
						flux0 = control.cFlux;
					}
					nMax++;
				}
				
				printf ("%lf %lf %lf %f %d\n", tdiff, fdiff, control.cFlux, acfStructure.probability, nMax);
				//fprintf (fin, "%lf %lf %lf %f\n", tdiff, fdiff, control.cFlux, acfStructure.probability);
			}
		}

		fflush (stdout);
		MPI_Finalize ();

		printf ("Threshold: %f \n", noiseStructure.detection);
		// deallocate memory
		deallocateMemory (&acfStructure, &noiseStructure);
	}
	//else
	//{
	//	plotDynSpec(pname, dname);
	//}

	/*
	if (fclose(fin))
	{
		printf ("Can't close file...\n");
		exit(1);
	}
	*/

	free(tdiss);
	free(fdiss);
	return 0;
}

