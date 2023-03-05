/*
File: phmmtools-common.hpp
Package: phmmtools
Copyright (C) 2008-2009 Aditya Bhargava (abhargava@cs.ualberta.ca)

phmmtools is free software: you can redistribute it and/or modify
it only for non-commercial use under the terms of the GNU General
Public License version 3 as published by the Free Software
Foundation.

phmmtools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with phmmtools.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PHMMTOOLS_COMMON_HPP
#define PHMMTOOLS_COMMON_HPP

#include <cassert>
#include <cmath>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#define ALPHABET 26

using namespace std;

unsigned int ctoi(char c);

unsigned int dtoi(double d);

char itoc(unsigned int i);

bool readClusters(string &file, vector< string > &labels, vector< vector< vector< unsigned int > > > &clusters);

bool readCharFreqs(string &file, double* q);

bool readCPs(string &file, double** cp);

bool writeAlignments(string &file, vector< string > &labels, vector< vector< vector< unsigned int > > > &alignments);

#endif
