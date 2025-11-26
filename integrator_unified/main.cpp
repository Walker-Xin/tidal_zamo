#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <time.h>
#include <complex>
#include "Eigen/Eigenvalues"

using namespace std;

#include "def.h"
#include "overloads.cpp"
#include "metric.cpp"
#include "christoffel.cpp"
#include "Cloc_newman.cpp"
#include "Cloc_spherical.cpp"
#include "diffeqs.cpp"
#include "conversion.cpp"
#include "eigen_compute.cpp"
#include "raytracing.cpp"

int main(int argc, char *argv[])
{
	// Variable declarations
	double error_tolerance = 1e-10; // C++ double has 15 decimal digits of precision
	int total_iterations;
	int eigenswitch;

	char filename_o[128];

	// Set default parameters
	// Note: Do not set spin = 0 or lam = 1 as they lead to instability; modify metric in these cases
	spin = 0.5;				  // spin parameter (0 < spin)
	charge = 0.5;			  // charge parameter
	lam = 0.5;				  // (cosine of) inclination angle (0 < lam <= 1)
	x = 2.701418260779751;	  // IBSO radius
	total_iterations = 10000; // a run usually takes less than 10k iterations
	scale = 1;				  // scale factor for lz
	eigenswitch = 1;		  // switch for eigenvalue computation

	// Set parameters from user input if provided
	if (argc > 1 && argc < 8)
	{
		printf("Usage: %s <spin> <charge> <lam> <x> <total_iterations> <scale> <eigenswitch>\n", argv[0]);
		return 1;
	}
	else if (argc == 8)
	{
		spin = atof(argv[1]);			  // spin parameter
		charge = atof(argv[2]);			  // charge parameter
		lam = atof(argv[3]);			  // inclination angle
		x = atof(argv[4]);				  // x parameter
		total_iterations = atoi(argv[5]); // number of iterations
		scale = atof(argv[6]);			  // scale factor for lz
		eigenswitch = atoi(argv[7]);	  // switch for eigenvalue computation
	}

	// Check if spin = 0
	if (spin == 0)
	{
		printf("Spin is zero, falling back to spherical spacetime.\n");
		spherical = 1;
		lam = 0; // set lam to zero (theta = pi/2) since it does not matter in spherical spacetime
	}

	// Get inclination angle in radians
	psi = acos(lam);

	// Get conserved quantities
	energy = 1.0;
	ep = energy; // alias for energy

	// Calculate IBSO lz
	lz = (-spin * spin + x * (x - charge * charge) - sqrt(x) * (x * x - 2 * x + spin * spin + charge * charge)) / (spin * (x - 1));
	if (spherical == 1)
	{
		lz = sqrt(x * x * (2 * x - charge * charge) / (x * x - 2 * x + charge * charge));
	}

	// Check if lz is negative
	if (lz < 0)
	{
		printf("Retrograde orbit\n");
		retrograde = 1;
	}
	else
	{
		printf("Prograde orbit\n");
		retrograde = 0;
	}

	// Calculate carter constant q and original carter constant k
	carter = lz * lz * (cos(psi) * cos(psi)) / (sin(psi) * sin(psi));
	if (spherical == 1)
	{
		carter = 0;
	}
	carter_origin = carter + (lz - spin * energy) * (lz - spin * energy);

	// Scale lz if necessary
	if (scale != 1)
	{
		lz *= scale;
	}

	// Print information
	printf("--- Parameters ---\n");
	printf("spin = %f, charge = %f, lam = %f, x = %f, total_iterations = %d, scale = %f\n", spin, charge, lam, x, total_iterations, scale);

	// Check if operating in naked singularity regime
	if (spin * spin + charge * charge > 1)
	{
		printf("Naked singularity regime\n");
	}
	else
	{
		printf("Black hole regime\n");
	}

	printf("--- Conserved quantities ---\n");
	printf("energy = %f, lz = %f,  carter_q = %f, carter_k = %f\n", energy, lz, carter, carter_origin);

	// If carter_origin is negative, abort
	if (carter_origin < 0)
	{
		printf("Carter_k is negative, aborting\n");
		return 1;
	}

	// Create output name
	sprintf(filename_o, "data/trace_a_%.2f_Q_%.2f_lambda_%.2f_%s.dat", spin, charge, lam, retrograde ? "ret" : "pro");

	// Start timer
	clock_t start, end;
	start = clock();

	printf("=== START ===\n");
	raytrace(filename_o, total_iterations, error_tolerance);

	// End timer
	end = clock();
	printf("Time taken by simulation: %f\n", (double)(clock() - start) / CLOCKS_PER_SEC);

	if (eigenswitch == 0)
	{
		printf("Eigenvalue computation disabled\n");
		end = clock();
		printf("=== END ===\n\n");
		return 0;
	}

	// Eigenvalue computation
	start = clock();
	eigen_compute(filename_o);
	end = clock();

	printf("Time taken by eigenvalue computation: %f\n", (double)(end - start) / CLOCKS_PER_SEC);
	printf("=== END ===\n\n");
	return double(end - start) / CLOCKS_PER_SEC;
}