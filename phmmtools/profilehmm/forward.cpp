/*
File: forward.cpp
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

void ProfileHMM::forward(vector< unsigned int > &seq, double **M, double **I, double **D, bool useq)
{
	double *q;

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
			M[k][i] = 0.0;

		// Insert states are numbered 0 to L
		for (unsigned int k = 0; k <= L; ++k)
			I[k][i] = 0.0;

		// Delete states are numbered 1 to L
		for (unsigned int k = 0; k <= L; ++k)
			D[k][i] = 0.0;
	}

	// initialisation
	M[0][0] = 1.0;

	// recursion
	// do the i = 0 case first, because in that case only the delete states count
	D[1][0] = a[M0][D1];

	for (unsigned int k = 2, Mk = M1 + 1, Ik = I0 + 2, Dk = D1 + 1; k <= L; ++k, ++Mk, ++Ik, ++Dk)
	{
		D[k][0] = D[k - 1][0] * a[Dk - 1][Dk];
	}

	// i is the matrix index, index is the sequence index
	for (unsigned int i = 1, index = 0; i <= seq.size(); ++i, ++index)
	{
		const unsigned int I1 = I0 + 1; // for convenience below

		// first do the I0 state (M0 is already done (initialisation) and D0 doesn't exist)
		I[0][i] = (e[I0][seq[index]] / q[seq[index]]) * (M[0][i - 1] * a[M0][I0] +
														 I[0][i - 1] * a[I0][I0]);

		// now do the k = 1 states (because the M1 and D1 states depend on D0, which doesn't exist)
		// since right now D0 is equivalent to the END state, this COULD be included in the loop with no
		// effect, but isn't to keep the code safe for future modification
		M[1][i] = (e[M1][seq[index]] / q[seq[index]]) * (M[0][i - 1] * a[M0][M1] +
														 I[0][i - 1] * a[I0][M1]);

		I[1][i] = (e[I1][seq[index]] / q[seq[index]]) * (M[1][i - 1] * a[M1][I1] + 
														 I[1][i - 1] * a[I1][I1] +
														 D[1][i - 1] * a[D1][I1]);

		D[1][i] = M[0][i] * a[M0][D1] +
				  I[0][i] * a[I0][D1];

		// now loop over the remaining k (2 to L)
		for (unsigned int k = 2, Mk = M1 + 1, Ik = I0 + 2, Dk = D1 + 1; k <= L; ++k, ++Mk, ++Ik, ++Dk)
		{
			M[k][i] = (e[Mk][seq[index]] / q[seq[index]]) * (M[k - 1][i - 1] * a[Mk - 1][Mk] +
															 I[k - 1][i - 1] * a[Ik - 1][Mk] +
															 D[k - 1][i - 1] * a[Dk - 1][Mk]);

			I[k][i] = (e[Ik][seq[index]] / q[seq[index]]) * (M[k][i - 1] * a[Mk][Ik] +
															 I[k][i - 1] * a[Ik][Ik] +
															 D[k][i - 1] * a[Dk][Ik]);

			D[k][i] = M[k - 1][i] * a[Mk - 1][Dk] +
					  I[k - 1][i] * a[Ik - 1][Dk] +
					  D[k - 1][i] * a[Dk - 1][Dk];
		}
	}

	// termination
	M[L + 1][seq.size() + 1] = M[L][seq.size()] * a[ML][END] +
							   I[L][seq.size()] * a[IL][END] +
							   D[L][seq.size()] * a[DL][END];

	if (!useq)
	{
		delete q;
	}
}

// forwardInLogSpace() is a replica of forward(), except that it is in log space.
// for clarity, I've refrained from replacing 0 with LOGZERO, and instead show the call to eln()
void ProfileHMM::forwardInLogSpace(vector< unsigned int > &seq, elnobj **M, elnobj **I, elnobj **D, bool useq)
{
	double *q;

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
			M[k][i] = eln(0.0);

		// Insert states are numbered 0 to L
		for (unsigned int k = 0; k <= L; ++k)
			I[k][i] = eln(0.0);

		// Delete states are numbered 1 to L
		for (unsigned int k = 0; k <= L; ++k)
			D[k][i] = eln(0.0);
	}

	// initialisation
	M[0][0] = eln(1.0);

	// recursion
	// do the i = 0 case first, because in that case only the delete states count
	D[1][0] = eln(a[M0][D1]);

	for (unsigned int k = 2, Mk = M1 + 1, Ik = I0 + 2, Dk = D1 + 1; k <= L; ++k, ++Mk, ++Ik, ++Dk)
	{
		D[k][0] = D[k - 1][0] * eln(a[Dk - 1][Dk]);
	}

	// i is the matrix index, index is the sequence index
	for (unsigned int i = 1, index = 0; i <= seq.size(); ++i, ++index)
	{
		const unsigned int I1 = I0 + 1; // for convenience below

		// first do the I0 state (M0 is already done (initialisation) and D0 doesn't exist)
		I[0][i] = eln(e[I0][seq[index]] / q[seq[index]]) * (M[0][i - 1] * eln(a[M0][I0]) + I[0][i - 1] * eln(a[I0][I0]));

		// now do the k = 1 states (because the M1 and D1 states depend on D0, which doesn't exist)
		// since right now D0 is equivalent to the END state, this COULD be included in the loop with no
		// effect, but isn't to keep the code safe for future modification
		M[1][i] = eln(e[M1][seq[index]] / q[seq[index]]) * (M[0][i - 1] * eln(a[M0][M1]) +
															I[0][i - 1] * eln(a[I0][M1]));

		I[1][i] = eln(e[I1][seq[index]] / q[seq[index]]) * (M[1][i - 1] * eln(a[M1][I1]) + 
															I[1][i - 1] * eln(a[I1][I1]) +
															D[1][i - 1] * eln(a[D1][I1]));

		D[1][i] = M[0][i] * eln(a[M0][D1]) +
				  I[0][i] * eln(a[I0][D1]);

		// now loop over the remaining k (2 to L)
		for (unsigned int k = 2, Mk = M1 + 1, Ik = I0 + 2, Dk = D1 + 1; k <= L; ++k, ++Mk, ++Ik, ++Dk)
		{
			M[k][i] = eln(e[Mk][seq[index]] / q[seq[index]]) * (M[k - 1][i - 1] * eln(a[Mk - 1][Mk]) +
																I[k - 1][i - 1] * eln(a[Ik - 1][Mk]) +
																D[k - 1][i - 1] * eln(a[Dk - 1][Mk]));

			I[k][i] = eln(e[Ik][seq[index]] / q[seq[index]]) * (M[k][i - 1] * eln(a[Mk][Ik]) +
																I[k][i - 1] * eln(a[Ik][Ik]) +
																D[k][i - 1] * eln(a[Dk][Ik]));

			D[k][i] = M[k - 1][i] * eln(a[Mk - 1][Dk]) +
					  I[k - 1][i] * eln(a[Ik - 1][Dk]) +
					  D[k - 1][i] * eln(a[Dk - 1][Dk]);
		}
	}

	// termination
	M[L + 1][seq.size() + 1] = M[L][seq.size()] * eln(a[ML][END]) +
							   I[L][seq.size()] * eln(a[IL][END]) +
							   D[L][seq.size()] * eln(a[DL][END]);

	if (!useq)
	{
		delete q;
	}
}
