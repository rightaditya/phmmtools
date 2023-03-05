/*
File: viterbi.cpp
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

vector< unsigned int > ProfileHMM::viterbi(vector< unsigned int > &seq, elnobj **M, elnobj **I, elnobj **D, bool useq)
{
	unsigned int **ptrM, **ptrI, **ptrD;
	elnobj temp1, temp2, temp3, themax;
	double *q;

	ptrM = new unsigned int *[L + 2];
	ptrI = new unsigned int *[L + 1];
	ptrD = new unsigned int *[L + 1];

	for (unsigned int i = 0; i < L + 2; ++i)
		ptrM[i] = new unsigned int[seq.size() + 2];

	for (unsigned int i = 0; i < L + 1; ++i)
	{
		ptrI[i] = new unsigned int[seq.size() + 2];

		ptrD[i] = new unsigned int[seq.size() + 2];
	}

	if (useq)
	{
		q = Settings::q;

		for (unsigned int i = 0; i < S; ++i)
		{
			if (q[i] == 0.0) // just a safeguard, there shouldn't be any q-values set to 0 anyway
			{
				q[i] = 1.0;
			}
		}
	}
	else
	{
		q = new double[S];

		for (unsigned int i = 0; i < S; ++i)
			q[i] = 1.0;
	}

	// clear the matrices
	for (unsigned int i = 0; i <= seq.size() + 1; ++i)
	{
		// Match states are numbered 0 (begin) to L + 1 (end)
		for (unsigned int k = 0; k <= L + 1; ++k)
		{
			M[k][i] = eln(0.0);
			ptrM[k][i] = 0;
		}

		// Insert states are numbered 0 to L
		for (unsigned int k = 0; k <= L; ++k)
		{
			I[k][i] = eln(0.0);
			ptrI[k][i] = 0;
		}

		// Delete states are numbered 1 to L
		for (unsigned int k = 0; k <= L; ++k)
		{
			D[k][i] = eln(0.0);
			ptrD[k][i] = 0;
		}
	}

	// initialisation
	M[0][0] = eln(1.0);

	// recursion
	// do the i = 0 case first, because in that case only the delete states count
	D[1][0] = eln(a[M0][D1]);
	ptrD[1][0] = M0;

	for (unsigned int k = 2, Mk = M1 + 1, Ik = I0 + 2, Dk = D1 + 1; k <= L; ++k, ++Mk, ++Ik, ++Dk)
	{
		D[k][0] = D[k - 1][0] * eln(a[Dk - 1][Dk]);
		ptrD[k][0] = Dk - 1;
	}

	// i is the matrix index, index is the sequence index
	for (unsigned int i = 1, index = 0; i <= seq.size(); ++i, ++index)
	{
		const unsigned int I1 = I0 + 1; // for convenience below

		// first do the I0 state (M0 is already done (initialisation) and D0 doesn't exist)
		temp1 = M[0][i - 1] * eln(a[M0][I0]);
		temp2 = I[0][i - 1] * eln(a[I0][I0]);
		themax = max(temp1, temp2);
		I[0][i] = eln(e[I0][seq[index]] / q[seq[index]]) * themax;
		if (themax == temp1)
			ptrI[0][i] = M0;
		else
			ptrI[0][i] = I0;

		// now do the j = 1 states (because the M1 and D1 states depend on D0, which doesn't exist)
		// since right now D0 is equivalent to the END state, this COULD be included in the loop with no
		// effect, but isn't to keep the code safe for future modification
		temp1 = M[0][i - 1] * eln(a[M0][M1]);
		temp2 = I[0][i - 1] * eln(a[I0][M1]);
		themax = max(temp1, temp2);
		M[1][i] = eln(e[M1][seq[index]] / q[seq[index]]) * themax;
		if (themax == temp1)
			ptrM[1][i] = M0;
		else
			ptrM[1][i] = I0;

		temp1 = M[1][i - 1] * eln(a[M1][I1]);
		temp2 = I[1][i - 1] * eln(a[I1][I1]);
		temp3 = D[1][i - 1] * eln(a[D1][I1]);
		themax = max(temp1, temp2, temp3);
		I[1][i] = eln(e[I1][seq[index]] / q[seq[index]]) * themax;
		if (themax == temp1)
			ptrI[1][i] = M1;
		else if (themax == temp2)
			ptrI[1][i] = I1;
		else
			ptrI[1][i] = D1;

		temp1 = M[0][i] * eln(a[M0][D1]);
		temp2 = I[0][i] * eln(a[I0][D1]);
		themax = max(temp1, temp2);
		D[1][i] = themax;
		if (themax == temp1)
			ptrD[1][i] = M0;
		else
			ptrD[1][i] = I0;

		// now loop over the remaining j (2 to L)
		for (unsigned int j = 2, Mj = M1 + 1, Ij = I0 + 2, Dj = D1 + 1; j <= L; ++j, ++Mj, ++Ij, ++Dj)
		{
			temp1 = M[j - 1][i - 1] * eln(a[Mj - 1][Mj]);
			temp2 = I[j - 1][i - 1] * eln(a[Ij - 1][Mj]);
			temp3 = D[j - 1][i - 1] * eln(a[Dj - 1][Mj]);
			themax = max(temp1, temp2, temp3);
			M[j][i] = eln(e[Mj][seq[index]] / q[seq[index]]) * themax;
			if (themax == temp1)
				ptrM[j][i] = Mj - 1;
			else if (themax == temp2)
				ptrM[j][i] = Ij - 1;
			else
				ptrM[j][i] = Dj - 1;

			temp1 = M[j][i - 1] * eln(a[Mj][Ij]);
			temp2 = I[j][i - 1] * eln(a[Ij][Ij]);
			temp3 = D[j][i - 1] * eln(a[Dj][Ij]);
			themax = max(temp1, temp2, temp3);
			I[j][i] = eln(e[Ij][seq[index]] / q[seq[index]]) * themax;
			if (themax == temp1)
				ptrI[j][i] = Mj;
			else if (themax == temp2)
				ptrI[j][i] = Ij;
			else
				ptrI[j][i] = Dj;

			temp1 = M[j - 1][i] * eln(a[Mj - 1][Dj]);
			temp2 = I[j - 1][i] * eln(a[Ij - 1][Dj]);
			temp3 = D[j - 1][i] * eln(a[Dj - 1][Dj]);
			themax = max(temp1, temp2, temp3);
			D[j][i] = themax;
			if (themax == temp1)
				ptrD[j][i] = Mj - 1;
			else if (themax == temp2)
				ptrD[j][i] = Ij - 1;
			else
				ptrD[j][i] = Dj - 1;
		}
	}

	// termination
	temp1 = M[L][seq.size()] * eln(a[ML][END]);
	temp2 = I[L][seq.size()] * eln(a[IL][END]);
	temp3 = D[L][seq.size()] * eln(a[DL][END]);
	themax = max(temp1, temp2, temp3);
	M[L + 1][seq.size() + 1] = themax;
	if (themax == temp1)
		ptrM[L + 1][seq.size() + 1] = ML;
	else if (themax == temp2)
		ptrM[L + 1][seq.size() + 1] = IL;
	else
		ptrM[L + 1][seq.size() + 1] = DL;

	// traceback
	vector< unsigned int > path;

	path.push_back(ptrM[L + 1][seq.size() + 1]);
	unsigned int i = seq.size();
	while (path.back() != M0)
	{
		unsigned int state = path.back();

		if (isM(state))
			path.push_back(ptrM[state - M0][i]);
		else if (isI(state))
			path.push_back(ptrI[state - I0][i]);
		else if (isD(state))
		{
			path.push_back(ptrD[state - D1 + 1][i]);

			++i; // because delete states are silent
		}

		--i;
	}

	path.pop_back(); // don't care about begin state

	reverse(path.begin(), path.end());

	// cleanup
	if (!useq)
	{
		delete q;
	}

	for (unsigned int i = 0; i < L + 2; ++i)
	{
		delete ptrM[i];
	}

	for (unsigned int i = 0; i < L + 1; ++i)
	{
		delete ptrI[i];
		delete ptrD[i];
	}

	delete ptrM; delete ptrI; delete ptrD;

	return path;
}

// The alignments returned can take on the following values:
//  -0 to S-1 to represent a character
//  -S to represent a deletion
//  -S + 1 to represent an insertion
//  -S + 2 is used in the guide state for a match column
//  -S + 3 is used in the guide state for an insert column
vector< vector< unsigned int > > ProfileHMM::multiAlign(vector< vector< unsigned int > > &seqsRaw)
{
	vector< vector< unsigned int > > seqs(seqsRaw.size());
	vector< vector< unsigned int > > paths(seqs.size());
	unsigned int *maxIVisits = new unsigned int[L + 1]; // for I0...IL

	for (unsigned int i = 0; i <= L; ++i)
		maxIVisits[i] = 0;

	for (unsigned int j = 0; j < seqsRaw.size(); ++j)
	{
		vector< unsigned int > seqRaw(seqsRaw[j]);
		vector< unsigned int > seq;

		// strip the input sequence of any inserts (i.e. clean it up)
		for (unsigned int i = 0; i < seqRaw.size(); ++i)
			if (seqRaw[i] < S)
				seq.push_back(seqRaw[i]);

		seqs[j] = seq;
	}

	for (unsigned int j = 0; j < seqs.size(); ++j)
	{
		elnobj **vM, **vI, **vD;
		vector< unsigned int > seq(seqs[j]);
		unsigned int *iVisits = new unsigned int[L + 1];

		vM = new elnobj *[L + 2];
		vI = new elnobj *[L + 1];
		vD = new elnobj *[L + 1];

		for (unsigned int i = 0; i < L + 2; ++i)
			vM[i] = new elnobj[seq.size() + 2];

		for (unsigned int i = 0; i < L + 1; ++i)
		{
			vI[i] = new elnobj[seq.size() + 2];

			vD[i] = new elnobj[seq.size() + 2];
		}

		for (unsigned int i = 0; i <= L; ++i)
			iVisits[i] = 0;		

		paths[j] = viterbi(seq, vM, vI, vD, true);

		for (unsigned int i = 0; i < paths[j].size(); ++i)
			if (isI(paths[j][i]))
				++iVisits[paths[j][i] - I0];

		for (unsigned int i = 0; i <= L; ++i)
			maxIVisits[i] = max(maxIVisits[i], iVisits[i]);

		// cleanup
		for (unsigned int i = 0; i < L + 2; ++i)
			delete vM[i];

		for (unsigned int i = 0; i < L + 1; ++i)
		{
			delete vI[i];
			delete vD[i];
		}

		delete vM; delete vI; delete vD;

		delete iVisits;
	}

	// now maxIVisits[i] contains the highest number of visits in a single path to state Ii

	vector< vector< unsigned int > > alignments(seqs.size() + 1);
	vector< unsigned int > guide;

	for (unsigned int j = 0; j < maxIVisits[0]; ++j)
		guide.push_back(S + 3);

	for (unsigned int i = 1; i <= L; ++i)
	{
		guide.push_back(S + 2);

		for (unsigned int j = 0; j < maxIVisits[i]; ++j)
			guide.push_back(S + 3);
	}
	alignments[0] = guide;

	for (unsigned int i = 0; i < paths.size(); ++i)
	{
		vector< unsigned int > path = paths[i];
		vector< unsigned int > seq = seqs[i];
		vector< unsigned int > alignment;
		vector< unsigned int > mBuffer(L + 1);
		vector< vector< unsigned int > > iBuffer(L + 1);

		for (unsigned int j = 0, k = 0; j < path.size(); ++j, ++k)
		{
			if (isM(path[j]))
				mBuffer[path[j] - M0] = seq[k];
			else if (isI(path[j]))
				iBuffer[path[j] - I0].push_back(seq[k]);
			else if (isD(path[j]))
			{
				mBuffer[path[j] - D1 + 1] = S;

				--k; // because delete states are silent
			}
		}

		for (unsigned int j = 0; j <= L; ++j)
		{
			int insertChars = maxIVisits[j] - iBuffer[j].size();

			for (int k = 0; k < insertChars; ++k)
				iBuffer[j].push_back(S + 1);
		}

		for (unsigned int k = 0; k < iBuffer[0].size(); ++k)
			alignment.push_back(iBuffer[0][k]);
		for (unsigned int j = 1; j <= L; ++j)
		{
			alignment.push_back(mBuffer[j]);

			for (unsigned int k = 0; k < iBuffer[j].size(); ++k)
				alignment.push_back(iBuffer[j][k]);
		}

		alignments[i + 1] = alignment;
	}

	delete maxIVisits;

	return alignments;
}
