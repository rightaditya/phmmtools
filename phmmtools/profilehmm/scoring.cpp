/*
File: scoring.cpp
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

double ProfileHMM::loglikelihood(vector< vector< unsigned int > > &seqsRaw)
{
	elnobj logP = eln(0.0);

	vector< vector< unsigned int > > seqs;
	for (unsigned int i = 0; i < seqsRaw.size(); ++i)
	{
		vector< unsigned int > seqRaw(seqsRaw[i]);
		vector< unsigned int > seq;

		// strip the input sequence of any inserts (i.e. clean it up)
		for (unsigned int j = 0; j < seqRaw.size(); ++j)
			if (seqRaw[j] < S)
				seq.push_back(seqRaw[j]);

		seqs.push_back(seq);
	}

	for (unsigned int j = 0; j < seqs.size(); ++j)
	{
		vector< unsigned int > seq(seqs[j]);
		elnobj **fM, **fI, **fD;

		fM = new elnobj *[L + 2];
		fI = new elnobj *[L + 1];
		fD = new elnobj *[L + 1];

		for (unsigned int i = 0; i < L + 2; ++i)
			fM[i] = new elnobj[seq.size() + 2];

		for (unsigned int i = 0; i < L + 1; ++i)
		{
			fI[i] = new elnobj[seq.size() + 2];
			fD[i] = new elnobj[seq.size() + 2];
		}

		forwardInLogSpace(seq, fM, fI, fD, false);

		logP += fM[L + 1][seq.size() + 1];

		// cleanup
		for (unsigned int i = 0; i < L + 2; ++i)
			delete fM[i];

		for (unsigned int i = 0; i < L + 1; ++i)
		{
			delete fI[i];
			delete fD[i];
		}

		delete fM; delete fI; delete fD;
	}

	return logP.value;
}

double ProfileHMM::logodds(vector< unsigned int > seqRaw)
{
	elnobj logP = eln(0.0);
	elnobj **fM, **fI, **fD;

	// strip the input sequence of any inserts (i.e. clean it up)
	vector< unsigned int > seq;
	for (unsigned int i = 0; i < seqRaw.size(); ++i)
		if (seqRaw[i] < S)
			seq.push_back(seqRaw[i]);

	fM = new elnobj *[L + 2];
	fI = new elnobj *[L + 1];
	fD = new elnobj *[L + 1];

	for (unsigned int i = 0; i < L + 2; ++i)
		fM[i] = new elnobj[seq.size() + 2];

	for (unsigned int i = 0; i < L + 1; ++i)
	{
		fI[i] = new elnobj[seq.size() + 2];
		fD[i] = new elnobj[seq.size() + 2];
	}

	forwardInLogSpace(seq, fM, fI, fD, true);

	logP = fM[L + 1][seq.size() + 1];

	// cleanup
	for (unsigned int i = 0; i < L + 2; ++i)
		delete fM[i];

	for (unsigned int i = 0; i < L + 1; ++i)
	{
		delete fI[i];
		delete fD[i];
	}

	delete fM; delete fI; delete fD;

	return logP.value;
}
