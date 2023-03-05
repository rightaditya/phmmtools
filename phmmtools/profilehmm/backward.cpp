/*
File: backward.cpp
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

void ProfileHMM::backward(vector< unsigned int > &seq, double **M, double **I, double **D, bool useq)
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
	M[L + 1][seq.size() + 1] = 1.0;
	M[L][seq.size()] = a[ML][END];
	I[L][seq.size()] = a[IL][END];
	D[L][seq.size()] = a[DL][END];

	// recursion
	// do the i = seq.size() case first, because in that case only the delete states count
	// this isn't included in the loop below to avoid the k = L case, which is really the initialisation done above
	// loop k from L - 1 down to 1
	for (int k = L - 1, Mk = ML - 1, Ik = IL - 1, Dk = DL - 1; k >= 1; --k, --Mk, --Ik, --Dk)
	{
		M[k][seq.size()] = D[k + 1][seq.size()] * a[Mk][Dk + 1];

		I[k][seq.size()] = D[k + 1][seq.size()] * a[Ik][Dk + 1];

		D[k][seq.size()] = D[k + 1][seq.size()] * a[Dk][Dk + 1];
	}

	// now the k = 0 case (not included in loop because D0 doesn't exist)
	M[0][seq.size()] = D[1][seq.size()] * a[M0][D1];
	I[0][seq.size()] = D[1][seq.size()] * a[I0][D1];

	// now loop i from seq.size() - 1 down to 0
	for (unsigned int i = seq.size() - 1; i >= 0; --i)
	{
		// k = L case (not included in loop because there is no D[L + 1])
		M[L][i] = I[L][i + 1] * a[ML][IL] * e[IL][seq[i]] / q[seq[i]];

		I[L][i] = I[L][i + 1] * a[IL][IL] * e[IL][seq[i]] / q[seq[i]];

		D[L][i] = I[L][i + 1] * a[DL][IL] * e[IL][seq[i]] / q[seq[i]];

		// loop k from L - 1 down to 1
		for (int k = L - 1, Mk = ML - 1, Ik = IL - 1, Dk = DL - 1; k >= 1; --k, --Mk, --Ik, --Dk)
		{
			M[k][i] = M[k + 1][i + 1] * a[Mk][Mk + 1] * e[Mk + 1][seq[i]] / q[seq[i]] +
					  I[k][i + 1] * a[Mk][Ik] * e[Ik][seq[i]] / q[seq[i]] +
					  D[k + 1][i] * a[Mk][Dk + 1];

			I[k][i] = M[k + 1][i + 1] * a[Ik][Mk + 1] * e[Mk + 1][seq[i]] / q[seq[i]] +
					  I[k][i + 1] * a[Ik][Ik] * e[Ik][seq[i]] / q[seq[i]] +
					  D[k + 1][i] * a[Ik][Dk + 1];

			D[k][i] = M[k + 1][i + 1] * a[Dk][Mk + 1] * e[Mk + 1][seq[i]] / q[seq[i]] +
					  I[k][i + 1] * a[Dk][Ik] * e[Ik][seq[i]] / q[seq[i]] +
					  D[k + 1][i] * a[Dk][Dk + 1];
		}

		// now the k = 0 case (not included in loop because D0 doesn't exist)
		M[0][i] = M[1][i + 1] * a[M0][M1] * e[M1][seq[i]] / q[seq[i]] +
				  I[0][i + 1] * a[M0][I0] * e[I0][seq[i]] / q[seq[i]] +
				  D[1][i] * a[M0][D1];

		I[0][i] = M[1][i + 1] * a[I0][M1] * e[M1][seq[i]] / q[seq[i]] +
				  I[0][i + 1] * a[I0][I0] * e[I0][seq[i]] / q[seq[i]] +
				  D[1][i] * a[I0][D1];
	}

	// termination
	M[0][0] = M[1][1] * a[M0][M1] * e[M1][seq[0]] / q[seq[0]] +
			  I[0][1] * a[M0][I0] * e[I0][seq[0]] / q[seq[0]] +
			  D[1][0] * a[M0][D1];

	if (!useq)
	{
		delete q;
	}
}

// backwardInLogSpace() is a replica of backward(), except that it is in log space.
// for clarity, I've refrained from replacing 0 with LOGZERO, and instead show the call to eln()
void ProfileHMM::backwardInLogSpace(vector< unsigned int > &seq, elnobj **M, elnobj **I, elnobj **D, bool useq)
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
	M[L + 1][seq.size() + 1] = eln(1.0);
	M[L][seq.size()] = eln(a[ML][END]);
	I[L][seq.size()] = eln(a[IL][END]);
	D[L][seq.size()] = eln(a[DL][END]);

	// recursion
	// do the i = seq.size() case first, because in that case only the delete states count
	// loop k from L - 1 down to 1
	for (int k = L - 1, Mk = ML - 1, Ik = IL - 1, Dk = DL - 1; k >= 1; --k, --Mk, --Ik, --Dk)
	{
		M[k][seq.size()] = D[k + 1][seq.size()] * eln(a[Mk][Dk + 1]);

		I[k][seq.size()] = D[k + 1][seq.size()] * eln(a[Ik][Dk + 1]);

		D[k][seq.size()] = D[k + 1][seq.size()] * eln(a[Dk][Dk + 1]);
	}

	// now the k = 0 case (not included in loop because D0 doesn't exist)
	M[0][seq.size()] = D[1][seq.size()] * eln(a[M0][D1]);
	I[0][seq.size()] = D[1][seq.size()] * eln(a[I0][D1]);

	// now loop i from seq.size() - 1 down to 0
	for (int i = seq.size() - 1; i >= 0; --i)
	{
		// k = L case (not included in loop because there is no D[L + 1])
		M[L][i] = I[L][i + 1] * eln(a[ML][IL]) * eln(e[IL][seq[i]] / q[seq[i]]);

		I[L][i] = I[L][i + 1] * eln(a[IL][IL]) * eln(e[IL][seq[i]] / q[seq[i]]);

		D[L][i] = I[L][i + 1] * eln(a[DL][IL]) * eln(e[IL][seq[i]] / q[seq[i]]);

		// loop k from L - 1 down to 1
		for (int k = L - 1, Mk = ML - 1, Ik = IL - 1, Dk = DL - 1; k >= 1; --k, --Mk, --Ik, --Dk)
		{
			M[k][i] = M[k + 1][i + 1] * eln(a[Mk][Mk + 1]) * eln(e[Mk + 1][seq[i]] / q[seq[i]]) +
					  I[k][i + 1] * eln(a[Mk][Ik]) * eln(e[Ik][seq[i]] / q[seq[i]]) +
					  D[k + 1][i] * eln(a[Mk][Dk + 1]);

			I[k][i] = M[k + 1][i + 1] * eln(a[Ik][Mk + 1]) * eln(e[Mk + 1][seq[i]] / q[seq[i]]) +
					  I[k][i + 1] * eln(a[Ik][Ik]) * eln(e[Ik][seq[i]] / q[seq[i]]) +
					  D[k + 1][i] * eln(a[Ik][Dk + 1]);

			D[k][i] = M[k + 1][i + 1] * eln(a[Dk][Mk + 1]) * eln(e[Mk + 1][seq[i]] / q[seq[i]]) +
					  I[k][i + 1] * eln(a[Dk][Ik]) * eln(e[Ik][seq[i]] / q[seq[i]]) +
					  D[k + 1][i] * eln(a[Dk][Dk + 1]);
		}

		// now the k = 0 case (not included in loop because D0 doesn't exist)
		M[0][i] = M[1][i + 1] * eln(a[M0][M1]) * eln(e[M1][seq[i]] / q[seq[i]]) +
				  I[0][i + 1] * eln(a[M0][I0]) * eln(e[I0][seq[i]] / q[seq[i]]) +
				  D[1][i] * eln(a[M0][D1]);

		I[0][i] = M[1][i + 1] * eln(a[I0][M1]) * eln(e[M1][seq[i]] / q[seq[i]]) +
				  I[0][i + 1] * eln(a[I0][I0]) * eln(e[I0][seq[i]] / q[seq[i]]) +
				  D[1][i] * eln(a[I0][D1]);
	}

	// termination
	M[0][0] = M[1][1] * eln(a[M0][M1]) * eln(e[M1][seq[0]] / q[seq[0]]) +
			  I[0][1] * eln(a[M0][I0]) * eln(e[I0][seq[0]] / q[seq[0]]) +
			  D[1][0] * eln(a[M0][D1]);

	if (!useq)
	{
		delete q;
	}
}
