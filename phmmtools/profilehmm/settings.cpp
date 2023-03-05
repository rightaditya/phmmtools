/*
File: settings.cpp
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

#include "settings.hpp"

// defaults, for allocation
unsigned long Settings::rngSeed = time(NULL);
double * Settings::q = NULL;
double ** Settings::cp = NULL;
pseudocount Settings::pcMethod = SUBSTITUTION;
double Settings::pcWeight = 5.0;

double Settings::initAlpha = 5.0;
bool Settings::initFavourM = true;
bool Settings::initInserts = true;
double Settings::initQScale = 3500;

double Settings::bwEpsilon = 1e-5;
unsigned int Settings::bwMaxIter = 0;
bool Settings::bwPCDuring = true;