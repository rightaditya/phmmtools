/*
File: state.cpp
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

vector< unsigned int > ProfileHMM::getDests(unsigned int src)
{
	vector< unsigned int > dsts;

	if (isM(src))
	{
		if (src == ML)
		{
			dsts.push_back(IL); // this (last) insert state
			dsts.push_back(END); // end state
		}
		else if (src < ML)
		{
			dsts.push_back(src - M0 + I0); // this insert state
			dsts.push_back(src - M0 + D1); // next delete state
			dsts.push_back(src + 1); // next match state
		}
	}
	else if (isI(src))
	{
		if (src == IL)
		{
			dsts.push_back(IL); // this (last) insert state
			dsts.push_back(END); // end state
		}
		else
		{
			dsts.push_back(src); // this insert state
			dsts.push_back(src - I0 + D1); // next delete state
			dsts.push_back(src - I0 + M1); // next match state
		}
	}
	else if (isD(src))
	{
		if (src == DL)
		{
			dsts.push_back(IL); // this (last) insert state
			dsts.push_back(END); // end state
		}
		else
		{
			dsts.push_back(src + 1); // next delete state
			dsts.push_back(src - D1 + I0 + 1); // this insert state
			dsts.push_back(src - D1 + M1 + 1); // next match state
		}
	}

	return dsts;
}

bool ProfileHMM::isM(unsigned int m)
{
	return (m >= M0) && (m <= END);
}

bool ProfileHMM::isI(unsigned int i)
{
	return (i >= I0) && (i <= IL);
}

bool ProfileHMM::isD(unsigned int d)
{
	return (d >= D1) && (d <= DL);
}
