/*
File: baum-welch.cpp
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

// Note: cp is a matrix used for substitution-matrix-based pseudocounts.  If it isn't set to NULL, q is ignored.
// It represents the conditional probabilities so that cp[i][j] = P(i|j).
int ProfileHMM::baumWelch(vector< vector< unsigned int> > &seqs)
{
	unsigned int iterations = 1;
	double oldlogP = loglikelihood(seqs), newlogP = baumWelchOnce(seqs, Settings::bwPCDuring);

	if (Settings::bwMaxIter != 0)
	{
		while (abs(newlogP - oldlogP) >= Settings::bwEpsilon && iterations < Settings::bwMaxIter)
		{
			oldlogP = newlogP;
			newlogP = baumWelchOnce(seqs, Settings::bwPCDuring);

			++iterations;
		}
	}
	else
	{
		while (abs(newlogP - oldlogP) >= Settings::bwEpsilon)
		{
			oldlogP = newlogP;
			newlogP = baumWelchOnce(seqs, Settings::bwPCDuring);

			++iterations;
		}
	}

	if (!Settings::bwPCDuring)
		newlogP = baumWelchOnce(seqs, true);

	return iterations;
}

// for clarity, I've refrained from replacing 0 with LOGZERO, and instead show the call to eln()
double ProfileHMM::baumWelchOnce(vector< vector< unsigned int> > &seqsRaw, bool addPseudocounts)
{
	elnobj **A, **E, *Asum, *Esum;

	// set up the A and E matrices
	A = new elnobj *[N];
	for (unsigned int i = 0; i < N; ++i)
		A[i] = new elnobj[N];
	E = new elnobj *[N];
	for (unsigned int i = 0; i < N; ++i)
		E[i] = new elnobj[S];
	Asum = new elnobj [N];
	Esum = new elnobj [N];

	// set all the A and E variables to zero
	for (unsigned int i = 0; i < N; ++i)
	{
		for (unsigned int j = 0; j < N; ++j)
			A[i][j] = eln(0.0);

		for (unsigned int j = 0; j < S; ++j)
			E[i][j] = eln(0.0);

		Asum[i] = eln(0.0);
		Esum[i] = eln(0.0);
	}

	// for each sequence j (i.e. for each sequence)
	for (unsigned int j = 0; j < seqsRaw.size(); ++j)
	{
		vector< unsigned int > seqRaw(seqsRaw[j]);
		vector< unsigned int > seq;

		for (unsigned int i = 0; i < seqRaw.size(); ++i)
			if (seqRaw[i] < S)
				seq.push_back(seqRaw[i]);

		elnobj **fM, **fI, **fD, **bM, **bI, **bD;

		fM = new elnobj *[L + 2]; bM = new elnobj *[L + 2];
		fI = new elnobj *[L + 1]; bI = new elnobj *[L + 1];
		fD = new elnobj *[L + 1]; bD = new elnobj *[L + 1];

		for (unsigned int i = 0; i < L + 2; ++i)
		{
			fM[i] = new elnobj[seq.size() + 2];
			bM[i] = new elnobj[seq.size() + 2];
		}

		for (unsigned int i = 0; i < L + 1; ++i)
		{
			fI[i] = new elnobj[seq.size() + 2];
			bI[i] = new elnobj[seq.size() + 2];

			fD[i] = new elnobj[seq.size() + 2];
			bD[i] = new elnobj[seq.size() + 2];
		}

		// do forward & backward algorithms
		forwardInLogSpace(seq, fM, fI, fD, false);
		backwardInLogSpace(seq, bM, bI, bD, false);

		elnobj logP = fM[L + 1][seq.size() + 1];

		// emission counts
		// for each symbol s in the output alphabet
		for (unsigned int s = 0; s < S; ++s)
		{
			// k = 0 case (since M0 has no emissions, we do I0 separately)
			elnobj bufI0 = eln(0.0);

			for (unsigned int i = 1; i <= seq.size(); ++i)
				if (seq[i - 1] == s)
					bufI0 += fI[0][i] * bI[0][i];

			bufI0 /= logP;
			E[I0][s] += bufI0;
			Esum[I0] += bufI0;

			// now loop from k = 1 to L
			for (unsigned int Mk = M1, Ik = I0 + 1, k = 1; Mk <= ML; ++Mk, ++Ik, ++k)
			{
				elnobj bufM = eln(0.0); elnobj bufI = eln(0.0);

				for (unsigned int i = 1; i <= seq.size(); ++i)
					if (seq[i - 1] == s)
					{
						bufM += fM[k][i] * bM[k][i];
						bufI += fI[k][i] * bI[k][i];
					}

				bufM /= logP; bufI /= logP;
				E[Mk][s] += bufM; E[Ik][s] += bufI;
				Esum[Mk] += bufM; Esum[Ik] += bufI;
			}
		}

		// transition counts
		for (unsigned int k = 0, n = 1, Mk = M0, Ik = I0, Dk = D1 - 1, Mn = M1;
			k <= L; ++k, ++n, ++Mk, ++Ik, ++Dk, ++Mn)
		{
			elnobj bufM = eln(0.0), bufI = eln(0.0), bufD = eln(0.0);
			bool DOK = Dk >= D1;
			unsigned int i;

			for (i = 0; i < seq.size(); ++i)
			{
				bufM += fM[k][i] * eln(a[Mk][Mn]) * eln(e[Mn][seq[i]]) * bM[n][i + 1];
				bufI += fI[k][i] * eln(a[Ik][Mn]) * eln(e[Mn][seq[i]]) * bM[n][i + 1];
				if (DOK)
					bufD += fD[k][i] * eln(a[Dk][Mn]) * eln(e[Mn][seq[i]]) * bM[n][i + 1];
			}

			// the i = seq.size() case is handled separately:
			// 1) because seq[seq.size()] doesn't exist
			// 2) note e[X][seq.size()] = 0 for all X except END, so we add the terms only in that case
			if (Mn == END)
			{
				bufM += fM[k][i] * eln(a[Mk][Mn]) * bM[n][i + 1];
				bufI += fI[k][i] * eln(a[Ik][Mn]) * bM[n][i + 1];
				if (DOK)
					bufD += fD[k][i] * eln(a[Dk][Mn]) * bM[n][i + 1];
			}

			bufM /= logP; bufI /= logP;
			A[Mk][Mn] += bufM; A[Ik][Mn] += bufI;
			Asum[Mk] += bufM; Asum[Ik] += bufI;

			if (DOK)
			{
				bufD /= logP;
				A[Dk][Mn] += bufD;
				Asum[Dk] += bufD;
			}
		}

		for (unsigned int k = 0, Mk = M0, Ik = I0, Dk = D1 - 1; k <= L; ++k, ++Mk, ++Ik, ++Dk)
		{
			elnobj bufM = eln(0.0), bufI = eln(0.0), bufD = eln(0.0);
			bool DOK = Dk >= D1;

			// because the i = seq.size() case gives zero results for anything but state END it can be
			// ignored here
			for (unsigned int i = 0; i < seq.size(); ++i)
			{
				bufM += fM[k][i] * eln(a[Mk][Ik]) * eln(e[Ik][seq[i]]) * bI[k][i + 1];
				bufI += fI[k][i] * eln(a[Ik][Ik]) * eln(e[Ik][seq[i]]) * bI[k][i + 1];
				if (DOK)
					bufD += fD[k][i] * eln(a[Dk][Ik]) * eln(e[Ik][seq[i]]) * bI[k][i + 1];
			}

			bufM /= logP; bufI /= logP;
			A[Mk][Ik] += bufM; A[Ik][Ik] += bufI;
			Asum[Mk] += bufM; Asum[Ik] += bufI;

			if (DOK)
			{
				bufD /= logP;
				A[Dk][Ik] += bufD;
				Asum[Dk] += bufD;
			}
		}

		for (unsigned int k = 0, n = 1, Mk = M0, Ik = I0, Dk = D1 - 1, Dn = D1;
			k < L; ++k, ++n, ++Mk, ++Ik, ++Dk, ++Dn)
		{
			elnobj bufM = eln(0.0), bufI = eln(0.0), bufD = eln(0.0);
			bool DOK = Dk >= D1;

			// because the i = seq.size() case gives zero results for anything but state END it can be
			// ignored here
			for (unsigned int i = 0; i < seq.size(); ++i)
			{
				bufM += fM[k][i] * eln(a[Mk][Dn]) * bD[n][i];
				bufI += fI[k][i] * eln(a[Ik][Dn]) * bD[n][i];
				if (DOK)
					bufD += fD[k][i] * eln(a[Dk][Dn]) * bD[n][i];
			}

			bufM /= logP; bufI /= logP;
			A[Mk][Dn] += bufM; A[Ik][Dn] += bufI;
			Asum[Mk] += bufM; Asum[Ik] += bufI;

			if (DOK)
			{
				bufD /= logP;
				A[Dk][Dn] += bufD;
				Asum[Dk] += bufD;
			}
		}

		// cleanup
		for (unsigned int i = 0; i < L + 2; ++i)
		{
			delete fM[i]; delete bM[i];
		}

		for (unsigned int i = 0; i < L + 1; ++i)
		{
			delete fI[i]; delete bI[i];
			delete fD[i]; delete bD[i];
		}

		delete fM; delete fI; delete fD;
		delete bM; delete bI; delete bD;
	}

	if (addPseudocounts)
		pseudocountifyInLogSpace(A, E, Asum, Esum);

	// build a
	for (unsigned int src = 0; src < N; ++src)
		if (eexp(Asum[src]) != 0)
			for (unsigned int dst = 0; dst < N; ++dst)
				a[src][dst] = eexp(A[src][dst] / Asum[src]);

	// build e
	for (unsigned int j = 0; j < S; ++j)
	{
		for (unsigned int Ik = I0; Ik <= IL; ++Ik)
			e[Ik][j] = eexp(E[Ik][j] / Esum[Ik]);

		for (unsigned int Mk = M1; Mk <= ML; ++Mk)
			e[Mk][j] = eexp(E[Mk][j] / Esum[Mk]);
	}

	// cleanup
	for (unsigned int i = 0; i < N; ++i)
	{
		delete A[i];
		delete E[i];
	}

	delete A; delete Asum;
	delete E; delete Esum;

	// return the new log likelihood of the model
	return loglikelihood(seqsRaw);
}
