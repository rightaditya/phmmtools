/*
File: phmmtools-common.cpp
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

#include "phmmtools-common.hpp"

unsigned int ctoi(char c)
{
	if (c == '-' || c == '.')
		return ALPHABET;
	else
		return toupper(c) - 0x41;
}

unsigned int dtoi(double d)
{
	return static_cast<int>(floor(d + 0.5));
}

char itoc(unsigned int i)
{
	if (i == ALPHABET)
		return '-'; // delete
	else if (i == ALPHABET + 1)
		return '-'; // insert
	else if (i == ALPHABET + 2)
		return 'M'; // match column for guide
	else if (i == ALPHABET + 3)
		return 'I'; // insert column for guide
	else
		return i + 0x41;
}

bool readClusters(string &file, vector< string > &labels, vector< vector< vector< unsigned int > > > &clusters)
{
	assert(clusters.empty()); assert(labels.empty()); // just in case
	ifstream fin(file.c_str());
	bool retVal = false, done = false;

	if (fin.is_open())
	{
		vector< vector< unsigned int > > cluster;

		while (!fin.eof() && !done)
		{
			string word;
			vector< unsigned int > wordInted;

			fin >> word;
			if (word.size() > 0 && word[0] != '#')
			{
				for (unsigned int i = 0; i < word.length(); ++i)
					wordInted.push_back(ctoi(word[i]));

				cluster.push_back(wordInted);
			}
			else if (word[0] == '#')
			{
				if (!cluster.empty())
				{
					clusters.push_back(cluster);

					cluster.clear();

					if (word[1] == '#')
						done = true;
				}

				if (!done)
					labels.push_back(word.substr(1));
			}
		}

		fin.close();
		
		retVal = true;
	}

	return retVal;
}

bool readCharFreqs(string &file, double* q)
{
	ifstream fin(file.c_str());
	bool retVal = false;

	if (fin.is_open())
	{
		unsigned int i = 0;

		while (!fin.eof() && i < ALPHABET)
		{
			fin >> q[i];
			++i;
		}

		fin.close();

		retVal = true;
	}

	return retVal;
}

bool readCPs(string &file, double** cp)
{
	ifstream fin(file.c_str());
	bool retVal = false;

	if (fin.is_open())
	{
		for (unsigned int i = 0; i < ALPHABET && !fin.eof(); ++i)
			for (unsigned int j = 0; j < ALPHABET && !fin.eof(); ++j)
				fin >> cp[i][j];

		fin.close();

		retVal = true;
	}

	return retVal;
}

bool writeAlignments(string &file, vector< string > &labels, vector< vector< vector< unsigned int > > > &alignments)
{
	assert(labels.size() == alignments.size()); // just in case
	ofstream fout(file.c_str(), ios::out | ios::trunc);
	bool retVal = false;

	if (fout.is_open())
	{
		for (unsigned int i = 0; i < labels.size(); ++i)
		{
			fout << "#" << labels[i] << endl;

			for (unsigned int j = 0; j < alignments[i].size(); ++j)
			{
				for (unsigned int k = 0; k < alignments[i][j].size(); ++k)
					fout << itoc(alignments[i][j][k]);
				
				fout << endl;
			}
		}

		fout << "##";

		retVal = true;
	}

	return retVal;
}
