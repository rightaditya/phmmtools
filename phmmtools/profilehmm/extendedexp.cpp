/*
File: extendedexp.cpp
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

#include "extendedexp.hpp"

elnobj operator + (const elnobj &lhs, const elnobj &rhs)
{
	elnobj ret;

	if (lhs.value == LOGZERO)
		ret.value = rhs.value;
	else if (rhs.value == LOGZERO)
		ret.value = lhs.value;
	else if (lhs.value > rhs.value)
		ret.value = lhs.value + eln(1.0 + exp(rhs.value - lhs.value)).value;
	else
		ret.value = rhs.value + eln(1.0 + exp(lhs.value - rhs.value)).value;

	return ret;
}

elnobj operator * (const elnobj &lhs, const elnobj &rhs)
{
	elnobj ret;

	if (lhs.value == LOGZERO || rhs.value == LOGZERO)
		ret.value = LOGZERO;
	else
		ret.value = lhs.value + rhs.value;

	return ret;
}

elnobj operator / (const elnobj &lhs, const elnobj &rhs)
{
	elnobj ret;
	
	// For now, just return LOGZERO when divisor is LOGZERO; not mathematically correct,
	// but good enough for use with profile HMM work.
	if (lhs.value == LOGZERO || rhs.value == LOGZERO)
		ret.value = LOGZERO;
	else
		ret.value = lhs.value - rhs.value;

	return ret;
}

void elnobj::operator += (const elnobj param)
{
	*this = *this + param;
}

void elnobj::operator *= (const elnobj param)
{
	*this = *this * param;
}

void elnobj::operator /= (const elnobj divisor)
{
	*this = *this / divisor;
}

bool operator != (const elnobj &lhs, const elnobj &rhs)
{
	return lhs.value != rhs.value;
}

bool operator == (const elnobj &lhs, const elnobj &rhs)
{
	return lhs.value == rhs.value;
}

bool operator > (const elnobj &lhs, const elnobj &rhs)
{
	return eexp(lhs) > eexp(rhs);
}

bool operator >= (const elnobj &lhs, const elnobj &rhs)
{
	return eexp(lhs) >= eexp(rhs);
}

bool operator < (const elnobj &lhs, const elnobj &rhs)
{
	return eexp(lhs) < eexp(rhs);
}

bool operator <= (const elnobj &lhs, const elnobj &rhs)
{
	return eexp(lhs) <= eexp(rhs);
}

double eexp(elnobj x)
{
	if (x.value == LOGZERO)
		return 0.0;
	else
		return exp(x.value);
}

elnobj eln(double x)
{
	elnobj ret;

	if (x == 0)
		ret.value = LOGZERO;
	else
		ret.value = log(x);

	return ret;
}

elnobj max(elnobj a, elnobj b, elnobj c)
{
	return max(a, max(b, c));
}
