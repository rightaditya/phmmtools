/*
File: maintenance.cpp
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

#include "profilehmm.hpp"

ProfileHMM::ProfileHMM(unsigned int length)
{
	L = length; N = 3 * L + 3; S = ALPHN;
	M0 = 0; M1 = 1; ML = L; END = L + 1;
	D1 = L + 2; DL = END + L; I0 = DL + 1; IL = I0 + L;

	rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng, Settings::rngSeed);

	initMatrices();

	initA();

	initE();
}

ProfileHMM::ProfileHMM(vector< vector< unsigned int > > alignedSeqs)
{
	bool allSeqsSameLength = true;

	for (unsigned int i = 0; i < alignedSeqs.size() - 1; ++i)
		if (alignedSeqs[i].size() != alignedSeqs[i + 1].size())
			allSeqsSameLength = false;

	assert(allSeqsSameLength); // just in case
	L = alignedSeqs[0].size(); S = ALPHN;
	vector< bool > matchState(alignedSeqs[0].size()), insertState(alignedSeqs[0].size());

	for (unsigned int i = 0; i < matchState.size(); ++i)
		matchState[i] = false;

	if (Settings::initInserts)
	{
		//MAPModelConstruct(alignedSeqs, matchState); // MAP doesn't work yet
		int thresh = (alignedSeqs.size() + 1) / 2;
		vector<int> counts(alignedSeqs[0].size());

		for (unsigned int i = 0; i < counts.size(); ++i)
			counts[i] = 0;

		for (unsigned int i = 0; i < alignedSeqs[0].size(); ++i)
			for (unsigned int j = 0; j < alignedSeqs.size(); ++j)
				if (alignedSeqs[j][i] < S)
					++counts[i];

		for (unsigned int i = 0; i < matchState.size(); ++i)
			matchState[i] = counts[i] >= thresh;
	}
	else
	{
		for (unsigned int i = 0; i < matchState.size(); ++i)
			matchState[i] = true;
	}

	for (unsigned int i = 0; i < matchState.size(); ++i)
		insertState[i] = !matchState[i];

	N = 3 * L + 3;
	M0 = 0; M1 = 1; ML = L; END = L + 1;
	D1 = L + 2; DL = END + L; I0 = DL + 1; IL = I0 + L;

	initMatrices();

	// build the general state sequence
	vector< unsigned int > stateSequence(alignedSeqs[0].size());
	
	// note the match states
	for (unsigned int i = 1, j = 0; i <= L && j < insertState.size(); ++i, ++j)
	{
		while (insertState[j])
			++j;

		stateSequence[j] = M0 + i;
	}

	// note the insert states
	if (insertState[0])
		stateSequence[0] = I0;

	for (unsigned int i = 1; i < insertState.size(); ++i)
	{
		if (insertState[i])
		{
			if (insertState[i - 1])
				stateSequence[i] = stateSequence[i - 1];
			else
				stateSequence[i] = stateSequence[i - 1] - M0 + I0;
		}
	}

	// now count & build
	double **A = new double *[N];
	double **E = new double *[N];
	double *Asum = new double[N];
	double *Esum = new double[N];
	
	for (unsigned int i = 0; i < N; ++i)
		A[i] = new double[N];
	for (unsigned int i = 0; i < N; ++i)
		E[i] = new double[S];

	// initialize As to 0
	for (unsigned int i = 0; i < N; ++i)
	{
		for (unsigned int j = 0; j < N; ++j)
			A[i][j] = 0.0;

		Asum[i] = 0.0;
	}

	// initialize Es to 0
	for (unsigned int i = 0; i < N; ++i)
	{
		for (unsigned int j = 0; j < S; ++j)
			E[i][j] = 0.0;

		Esum[i] = 0.0;
	}

	// count As
	for (unsigned int i = 0; i < alignedSeqs.size(); ++i)
	{
		vector< unsigned int > seqStateSequence;

		for (unsigned int j = 0; j < alignedSeqs[i].size(); ++j)
		{
			if (alignedSeqs[i][j] == S)
			{
				if (!insertState[j])
					seqStateSequence.push_back(stateSequence[j] - M1 + D1);
			}
			else
				seqStateSequence.push_back(stateSequence[j]);
		}

		unsigned int k = 0;

		++A[M0][seqStateSequence[k]];
		++Asum[M0];

		for (k = 0; k < seqStateSequence.size() - 1; ++k)
		{
			++A[seqStateSequence[k]][seqStateSequence[k + 1]];
			++Asum[seqStateSequence[k]];
		}

		++A[seqStateSequence[k]][END];
		++Asum[seqStateSequence[k]];
	}

	// count Es
	for (unsigned int i = 0; i < alignedSeqs.size(); ++i)
	{
		for (unsigned int j = 0; j < alignedSeqs[i].size(); ++j)
		{
			if (alignedSeqs[i][j] < S)
			{
				++E[stateSequence[j]][alignedSeqs[i][j]];
				++Esum[stateSequence[j]];
			}
		}
	}

	pseudocountify(A, E, Asum, Esum);

	// build As
	for (unsigned int i = 0; i < N; ++i)
		if (Asum[i] != 0.0)
			for (unsigned int j = 0; j < N; ++j)
				a[i][j] = A[i][j] / Asum[i];

	// build Es
	for (unsigned int j = 0; j < S; ++j)
	{
		for (unsigned int Ik = I0; Ik <= IL; ++Ik)
			if (Esum[Ik] != 0.0)
				e[Ik][j] = E[Ik][j] / Esum[Ik];

		for (unsigned int Mk = M1; Mk <= ML; ++Mk)
			if (Esum[Mk] != 0.0)
				e[Mk][j] = E[Mk][j] / Esum[Mk];
	}

	for (unsigned int i = 0; i < N; ++i)
	{
		delete A[i];
		delete E[i];
	}

	delete A; delete Asum;
	delete E; delete Esum;
}

ProfileHMM::~ProfileHMM()
{
	deleteMatrices();

	gsl_rng_free(rng);
}

void ProfileHMM::deleteMatrices()
{
	for (unsigned int i = 0; i < N; ++i)
	{
		delete a[i];
		delete e[i];
	}

	delete a;
	delete e;
}

// State divisions are defined as follows given length L:
// 0: Begin state
// 1 to L: Match states (L)
// L + 1: End state
// L + 2 to 2L + 1: Delete states (L)
// 2L + 2 to 3L + 2: Insert states (L + 1)
void ProfileHMM::initA()
{
	double alphas[3], probs[3];

	for (int i = 0; i < 3; ++i)
		alphas[i] = Settings::initAlpha;
	
	for (unsigned int Mi = M0, Mn = M1, Ii = I0, Di = D1 - 1, Dn = D1;
		 Mi < ML; ++Mi, ++Mn, ++Ii, ++Di, ++Dn)
	{
		// Match-outs
		gsl_ran_dirichlet(rng, 3, alphas, probs);
		if (Settings::initFavourM)
			swap(*max_element(probs, probs + 3), probs[0]);
		a[Mi][Mn] = probs[0];
		a[Mi][Ii] = probs[1];
		a[Mi][Dn] = probs[2];

		// Insert-outs
		gsl_ran_dirichlet(rng, 3, alphas, probs);
		if (Settings::initFavourM)
			swap(*max_element(probs, probs + 3), probs[0]);
		a[Ii][Mn] = probs[0];
		a[Ii][Ii] = probs[1];
		a[Ii][Dn] = probs[2];

		// Delete-outs (except D0 which doesn't exist)
		if (Di >= D1)
		{
			gsl_ran_dirichlet(rng, 3, alphas, probs);
			if (Settings::initFavourM)
				swap(*max_element(probs, probs + 3), probs[0]);
			a[Di][Mn] = probs[0];
			a[Di][Ii] = probs[1];
			a[Di][Dn] = probs[2];
		}
	}
	// since the last states of each type only have two transitions (instead of three), they are handled
	// separately
	gsl_ran_dirichlet(rng, 2, alphas, probs);
	if (Settings::initFavourM)
		swap(*max_element(probs, probs + 2), probs[0]);
	a[ML][END] = probs[0];
	a[ML][IL] = probs[1];

	gsl_ran_dirichlet(rng, 2, alphas, probs);
	if (Settings::initFavourM)
		swap(*max_element(probs, probs + 2), probs[0]);
	a[IL][END] = probs[0];
	a[IL][IL] = probs[1];

	gsl_ran_dirichlet(rng, 2, alphas, probs);
	if (Settings::initFavourM)
		swap(*max_element(probs, probs + 2), probs[0]);
	a[DL][END] = probs[0];
	a[DL][IL] = probs[1];
}

void ProfileHMM::initE()
{
	// set all insert-state emissions to the background frequencies
	for (unsigned int i = I0; i <= IL; ++i)
	{
		for (unsigned int j = 0; j < S; ++j)
			e[i][j] = Settings::q[j];
	}

	// sample match-state emissions from a Dirichlet distribution
	double* qscaled = new double[S];

	for (unsigned int i = 0; i < S; ++i)
		qscaled[i] = Settings::initQScale * Settings::q[i];

	for (unsigned int i = M1; i <= ML; ++i)
	{
		gsl_ran_dirichlet(rng, S, qscaled, e[i]);
	}

	delete qscaled;
}

void ProfileHMM::initMatrices()
{
	a = new double *[N];
	for (unsigned int i = 0; i < N; ++i)
		a[i] = new double[N];
	e = new double *[N];
	for (unsigned int i = 0; i < N; ++i)
		e[i] = new double[S];

	// initialize a to 0
	for (unsigned int i = 0; i < N; ++i)
		for (unsigned int j = 0; j < N; ++j)
			a[i][j] = 0.0;

	// initialize e to 1 (silent states are given emission probabilities of 1 for all output symbols)
	for (unsigned int i = 0; i < N; ++i)
		for (unsigned int j = 0; j < S; ++j)
			e[i][j] = 1.0;
}

// this doesn't work yet
void ProfileHMM::MAPModelConstruct(vector< vector< unsigned int > > alignedSeqs, vector< bool > matchState)
{
	// note that there is no penalty applied for longer models
	unsigned int mapL = alignedSeqs[0].size(); unsigned int arraySize = mapL + 2;
	double *score = new double[arraySize];
	unsigned int *sigma = new unsigned int[arraySize];
	double **T = new double*[arraySize];
	double *M = new double[arraySize];
	double **I = new double*[arraySize];
	for (unsigned int i = 0; i < arraySize; ++i)
	{
		T[i] = new double[arraySize];
		I[i] = new double[arraySize];

		score[i] = -1.0 * numeric_limits<double>::infinity();
		sigma[i] = 0;
		M[i] = 0.0;

		for (unsigned int j = 0; j < arraySize; ++j)
		{
			T[i][j] = 0.0;
			I[i][j] = 0.0;
		}
	}

	// variable setup
	// setup M - note no pseudocounts are used
	double *counts = new double[S], countsum;
	for (unsigned int i = 1; i <= mapL; ++i)
	{
		for (unsigned int j = 0; j < S; ++j)
			counts[j] = 0.0;
		countsum = 0.0;

		for (unsigned int j = 0; j < alignedSeqs.size(); ++j)
		{
			unsigned int val = alignedSeqs[j][i - 1];

			if (val < S)
			{
				++counts[val];
				++countsum;
			}
		}

		for (unsigned int j = 0; j < S; ++j)
			if (counts[j] > 0)
				M[i] += counts[j] * log(counts[j] / countsum);
	}

	// setup I - again, no pseudocounts
	for (unsigned int j = 0; j < arraySize; ++j)
	{
		for (unsigned int i = 0; i + 1 < j; ++i)
		{
			for (unsigned int k = 0; k < S; ++k)
				counts[k] = 0.0;
			countsum = 0.0;

			for (unsigned int k = 0; k < alignedSeqs.size(); ++k)
			{
				for (unsigned int l = i + 1; l < j; ++l)
				{
					unsigned int val = alignedSeqs[k][l - 1];

					if (val < S)
					{
						++counts[val];
						++countsum;
					}
				}
			}

			for (unsigned k = 0; k < S; ++k)
				if (counts[k] > 0)
					I[i + 1][j - 1] += counts[k] * log(counts[k] / countsum);

		}
	}
	delete counts;

	// setup T - again, no pseudocounts
	double transitions[3][3], transitionsum[3]; // indexed as 0: M, 1: I, 2: D
	for (unsigned int j = 1; j < arraySize; ++j)
	{
		for (unsigned int i = 0; i < j; ++i)
		{
			for (unsigned int k = 0; k < 3; ++k)
			{
				for (unsigned int l = 0; l < 3; ++l)
					transitions[k][l] = 0.0;
				transitionsum[k] = 0.0;
			}

			for (unsigned int k = 0; k < alignedSeqs.size(); ++k)
			{
				int numInserts = 0;
				bool srcM, dstM; // whether the "source" (i) and "destination" (j) columns for this sequence are match states or delete states

				for (unsigned int l = i + 1; l < j; ++l)
					if (alignedSeqs[k][l - 1] < S)
						++numInserts;

				if (i == 0)
					srcM = true;
				else if (alignedSeqs[k][i - 1] < S)
					srcM = true;
				else
					srcM = false;

				if (j == mapL + 1)
					dstM = true;
				else if (alignedSeqs[k][j - 1] < S)
					dstM = true;
				else
					dstM = false;

				if (numInserts > 0)
				{
					++transitions[srcM ? 0 : 2][1]; ++transitionsum[srcM ? 0 : 2]; // "source" to first insert

					transitions[1][1] += numInserts - 1; transitionsum[1] += numInserts - 1;

					++transitions[1][dstM ? 0 : 2]; ++transitionsum[1]; // last insert to "destination"
				}
				else
				{
					++transitions[srcM ? 0 : 2][dstM ? 0 : 2]; ++transitionsum[srcM ? 0 : 2]; // only one transition if there are no inserts
				}
			}

			for (unsigned int k = 0; k < 3; ++k)
				for (unsigned int l = 0; l < 3; ++l)
					if (transitions[k][l] > 0)
						T[i][j] += transitions[k][l] * log(transitions[k][l] / transitionsum[k]);
		}
	}

	// initialisation
	score[0] = 0.0;

	// recurrence
	for (unsigned int j = 1; j < arraySize; ++j)
	{
		for (unsigned int i = 0; i < j; ++i)
		{
			double trans = T[i][j];
			double mat = M[i];
			double ins = I[i + 1][j - 1];
			double sum = score[i] + trans + mat + ins;

			if (sum > score[j])
			{
				score[j] = sum;
				sigma[j] = i;
			}
		}
	}

	// traceback
	unsigned int j = sigma[mapL + 1];
	while (j > 0)
	{
		matchState[j - 1] = true;

		j = sigma[j];
	}

	for (unsigned int i = 0; i < arraySize; ++i)
	{
		delete I[i];
		delete T[i];
	}
	delete I;
	delete M;
	delete T;
	delete sigma;
	delete score;
}

void ProfileHMM::pseudocountify(double **A, double **E, double *Asum, double *Esum)
{
	// use Laplace's rule for transition pseudocounts
	for (unsigned int i = 0; i < N; ++i)
	{
		vector< unsigned int > dests = getDests(i);

		for (unsigned int j = 0; j < dests.size(); ++j)
		{
			A[i][dests[j]] += 1.0;
			Asum[i] += 1.0;
		}
	}

	if (Settings::pcMethod == CONSTANT)
	{
		for (unsigned int Ik = I0; Ik <= IL; ++Ik)
		{
			for (unsigned int s = 0; s < S; ++s)
				E[Ik][s] += Settings::pcWeight / static_cast<double>(S);

			Esum[Ik] += Settings::pcWeight;
		}
		
		for (unsigned int Mk = M1; Mk <= ML; ++Mk)
		{
			for (unsigned int s = 0; s < S; ++s)
				E[Mk][s] += Settings::pcWeight / static_cast<double>(S);

			Esum[Mk] += Settings::pcWeight;
		}
	}
	else if (Settings::pcMethod == BACKGROUND)
	{
		for (unsigned int Ik = I0; Ik <= IL; ++Ik)
		{
			for (unsigned int s = 0; s < S; ++s)
				E[Ik][s] += Settings::pcWeight * Settings::q[s];

			Esum[Ik] += Settings::pcWeight;
		}
		
		for (unsigned int Mk = M1; Mk <= ML; ++Mk)
		{
			for (unsigned int s = 0; s < S; ++s)
				E[Mk][s] += Settings::pcWeight * Settings::q[s];

			Esum[Mk] += Settings::pcWeight;
		}
	}
	else if (Settings::pcMethod == SUBSTITUTION)
	{
		double** alpha = new double*[N];

		for (unsigned int i = 0; i < N; ++i)
		{
			alpha[i] = new double[S];

			for (unsigned int j = 0; j < S; ++j)
				alpha[i][j] = 0.0;
		}

		for (unsigned int a = 0; a < S; ++a)
		{
			for (unsigned int Ik = I0; Ik <= IL; ++Ik)
			{
				if (Esum[Ik] != 0.0)
				{
					for (unsigned int b = 0; b < S; ++b)
					{
						alpha[Ik][a] += E[Ik][b] * Settings::cp[a][b] / Esum[Ik];
					}

					alpha[Ik][a] *= Settings::pcWeight;
				}
			}

			
			for (unsigned int Mk = M1; Mk <= ML; ++Mk)
			{
				if (Esum[Mk] != 0.0)
				{
					for (unsigned int b = 0; b < S; ++b)
					{
						alpha[Mk][a] += E[Mk][b] * Settings::cp[a][b] / Esum[Mk];
					}

					alpha[Mk][a] *= Settings::pcWeight;
				}
			}
		}

		for (unsigned int s = 0; s < S; ++s)
		{
			for (unsigned int Ik = I0; Ik <= IL; ++Ik)
			{
				E[Ik][s] += alpha[Ik][s];
				Esum[Ik] += alpha[Ik][s];
			}

			for (unsigned int Mk = M1; Mk <= ML; ++Mk)
			{
				E[Mk][s] += alpha[Mk][s];
				Esum[Mk] += alpha[Mk][s];
			}
		}

		for (unsigned int i = 0; i < N; ++i)
			delete alpha[i];

		delete alpha;
	}
}

// pseudocountifyInLogSpace() is a replica of pseudocountify(), except that it is in log space
void ProfileHMM::pseudocountifyInLogSpace(elnobj **A, elnobj **E, elnobj *Asum, elnobj *Esum)
{
	// use Laplace's rule for transition pseudocounts
	for (unsigned int i = 0; i < N; ++i)
	{
		vector< unsigned int > dests = getDests(i);

		for (unsigned int j = 0; j < dests.size(); ++j)
		{
			A[i][dests[j]] += eln(1.0);
			Asum[i] += eln(1.0);
		}
	}

	elnobj lnw = eln(Settings::pcWeight);

	if (Settings::pcMethod == CONSTANT)
	{
		elnobj lnS = eln(S);

		for (unsigned int Ik = I0; Ik <= IL; ++Ik)
		{
			for (unsigned int s = 0; s < S; ++s)
				E[Ik][s] += lnw / lnS;

			Esum[Ik] += lnw;
		}
		
		for (unsigned int Mk = M1; Mk <= ML; ++Mk)
		{
			for (unsigned int s = 0; s < S; ++s)
				E[Mk][s] += lnw / lnS;

			Esum[Mk] += lnw;
		}
	}
	else if (Settings::pcMethod == BACKGROUND)
	{
		elnobj *lnq = new elnobj[S];

		for (unsigned int i = 0; i < S; ++i)
			lnq[i] = eln(Settings::q[i]);

		for (unsigned int Ik = I0; Ik <= IL; ++Ik)
		{
			for (unsigned int s = 0; s < S; ++s)
				E[Ik][s] += lnw * lnq[s];

			Esum[Ik] += lnw;
		}
		
		for (unsigned int Mk = M1; Mk <= ML; ++Mk)
		{
			for (unsigned int s = 0; s < S; ++s)
				E[Mk][s] += lnw * lnq[s];

			Esum[Mk] += lnw;
		}

		delete lnq;
	}
	else if (Settings::pcMethod == SUBSTITUTION)
	{
		elnobj** alpha = new elnobj*[N];
		elnobj** lncp = new elnobj*[S];

		for (unsigned int i = 0; i < N; ++i)
		{
			alpha[i] = new elnobj[S];

			for (unsigned int j = 0; j < S; ++j)
				alpha[i][j] = eln(0.0);
		}

		for (unsigned int i = 0; i < S; ++i)
		{
			lncp[i] = new elnobj[S];

			for (unsigned int j = 0; j < S; ++j)
				lncp[i][j] = eln(Settings::cp[i][j]);
		}

		for (unsigned int a = 0; a < S; ++a)
		{
			for (unsigned int Ik = I0; Ik <= IL; ++Ik)
			{
				if (Esum[Ik] != eln(0.0))
				{
					for (unsigned int b = 0; b < S; ++b)
					{
						alpha[Ik][a] += E[Ik][b] * lncp[a][b] / Esum[Ik];
					}

					alpha[Ik][a] *= lnw;
				}
			}

			
			for (unsigned int Mk = M1; Mk <= ML; ++Mk)
			{
				if (Esum[Mk] != eln(0.0))
				{
					for (unsigned int b = 0; b < S; ++b)
					{
						alpha[Mk][a] += E[Mk][b] * lncp[a][b] / Esum[Mk];
					}

					alpha[Mk][a] *= lnw;
				}
			}
		}

		for (unsigned int s = 0; s < S; ++s)
		{
			for (unsigned int Ik = I0; Ik <= IL; ++Ik)
			{
				E[Ik][s] += alpha[Ik][s];
				Esum[Ik] += alpha[Ik][s];
			}

			for (unsigned int Mk = M1; Mk <= ML; ++Mk)
			{
				E[Mk][s] += alpha[Mk][s];
				Esum[Mk] += alpha[Mk][s];
			}
		}

		for (unsigned int i = 0; i < N; ++i)
			delete alpha[i];

		delete alpha;
	}
}
