/*
File: settings.hpp
Package: libprofilehmm
Copyright (C) 2008-2009 Aditya Bhargava (abhargava@cs.ualberta.ca)

libprofilehmm is free software: you can redistribute it and/or modify
it only for non-commercial use under the terms of the GNU General
Public License version 3 as published by the Free Software
Foundation.

libprofilehmm is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with libprofilehmm.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include <ctime>

using namespace std;

// alphabet size
// TODO: throw this in the main settings class at some point and make
// it properly parameterized
#define ALPHN 26

enum pseudocount
{
	CONSTANT, // a constant value added to each count (pseudocount)
	BACKGROUND, // a value proportional to the background distribution
	SUBSTITUTION // a substitution matrix-based method
};

// there are no checks whether the settings have been set or not before
// any profile HMM stuff is done, so that is a potential danger (TODO?)
struct Settings
{
	static unsigned long rngSeed; // seed for random number generator
	static double *q; // background frequencies
	static double **cp; // cp[a][b] = P[a|b]
	static pseudocount pcMethod; // pseudocount method
	static double pcWeight; // pseudocount weight

	// building a profile HMM from scratch
	static double initAlpha; // Dirichlet parameters for transitions
	static bool initFavourM; // whether or not to favour match states
	static bool initInserts; // whether or not to use insert states
	static double initQScale; // scaling for bg frequencies for Es

	// Baum-Welch settings
	static double bwEpsilon; // maximum difference between iterations
	static unsigned int bwMaxIter; // maximum number of iterations
	static bool bwPCDuring; // add pseudocounts between iterations
};

#endif
