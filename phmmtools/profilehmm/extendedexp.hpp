/*
File: extendedexp.hpp
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

#ifndef EXTENDEDEXP_HPP
#define EXTENDEDEXP_HPP

#include <algorithm>
#include <cmath>
#include <limits>

using namespace std;

#define LOGZERO (-1.0 * numeric_limits<double>::infinity())

// Note: operators overloaded only as needed
struct elnobj
{
	double value;

	void operator += (const elnobj param);
	void operator *= (const elnobj param);
	void operator /= (const elnobj divisor);
};

elnobj operator + (const elnobj &lhs, const elnobj &rhs);
elnobj operator * (const elnobj &lhs, const elnobj &rhs);
elnobj operator / (const elnobj &lhs, const elnobj &rhs);

bool operator == (const elnobj &lhs, const elnobj &rhs);
bool operator != (const elnobj &lhs, const elnobj &rhs);
bool operator > (const elnobj &lhs, const elnobj &rhs);
bool operator >= (const elnobj &lhs, const elnobj &rhs);
bool operator < (const elnobj &lhs, const elnobj &rhs);
bool operator <= (const elnobj &lhs, const elnobj &rhs);

double eexp(elnobj x);

elnobj eln(double x);

// For convenience of notation only
elnobj max(elnobj a, elnobj b, elnobj c);

#endif
