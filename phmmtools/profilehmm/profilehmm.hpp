/*
File: profilehmm.hpp
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

#ifndef PROFILEHMM_HPP
#define PROFILEHMM_HPP

#include <cassert>
#include <cmath>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "extendedexp.hpp"
#include "settings.hpp"

#ifndef NOMINMAX // from WinDef.h
#undef max
#endif

using namespace std;

class ProfileHMM
{
private:
	unsigned int L; // length of profile HMM, i.e., number of match states
	unsigned int N; // TOTAL number of states = L * 3 + 3
	unsigned int S; // size of output alphabet

	// probability tables
	double **a; // a[i][j] is the transition probability of going from state i to state j
	double **e; // e[i][k] is the probability of emitting symbol k while in state i

	// State divisions are defined as follows given length L:
	// 0: Begin state (labelled M0)
	unsigned int M0;
	// 1 to L: Match states (L)
	unsigned int M1, ML;
	// L + 1: End state
	unsigned int END;
	// L + 2 to 2L + 1: Delete states (L)
	unsigned int D1, DL;
	// 2L + 2 to 3L + 2: Insert states (L + 1)
	unsigned int I0, IL;

	// random number generator
	gsl_rng *rng;

	// init and delete stuff
	void deleteMatrices();
	void initA();
	void initE();
	void initMatrices();
	void MAPModelConstruct(vector< vector< unsigned int > > alignedSeqs, vector< bool > insertState);

	// state stuff
	vector< unsigned int > getDests(unsigned int state);
	bool isM(unsigned int m);
	bool isI(unsigned int i);
	bool isD(unsigned int d);

	double baumWelchOnce(vector< vector< unsigned int > > &seqsRaw, bool addPseudocounts);
	void pseudocountify(double **A, double **E, double *Asum, double *Esum);
	void pseudocountifyInLogSpace(elnobj **A, elnobj **E, elnobj *Asum, elnobj *Esum);

public:
	// constructor that builds new profile HMM with random parameters given length
	ProfileHMM(unsigned int length);
	ProfileHMM(vector< vector< unsigned int > > alignedSeqs);
	~ProfileHMM(); // destructor

	void backward(vector< unsigned int > &seq, double **M, double **I, double **D, bool useq);
	void backwardInLogSpace(vector< unsigned int > &seq, elnobj **M, elnobj **I, elnobj **D, bool useq);
	int baumWelch(vector< vector< unsigned int > > &seqs);
	void forward(vector< unsigned int > &seq, double **M, double **I, double **D, bool useq);
	void forwardInLogSpace(vector< unsigned int > &seq, elnobj **M, elnobj **I, elnobj **D, bool useq);
	double loglikelihood(vector< vector< unsigned int > > &seqsRaw);
	double logodds(vector< unsigned int > seqRaw);
	vector< vector< unsigned int > > multiAlign(vector< vector< unsigned int > > &seqsRaw);
	vector< unsigned int > viterbi(vector< unsigned int > &seq, elnobj **M, elnobj **I, elnobj **D, bool useq);
	void writeme(ostream& out);
};

#endif
