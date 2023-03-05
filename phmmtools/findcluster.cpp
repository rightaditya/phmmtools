/*
File: findcluster.cpp
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

#include <cassert>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "profilehmm/profilehmm.hpp"
#include "phmmtools-common.hpp"

using namespace std;

namespace po = boost::program_options;

void printScores(vector< string > labels, vector< double > scores);

int main(int argc, char *argv[])
{
	int retVal = 0;

	bool singles, prealigned;
	string clustersFile, qFile, word = "", cpFile;
	unsigned int pcMethod;

	vector< vector< vector< unsigned int > > > clusters;
	vector< string > labels;
	vector< unsigned int > wordInted;

	try
	{
		po::options_description generic("Generic options");
		generic.add_options()
			("help,h", "Prints help message.")
			;

		po::options_description config("Configuration options");
		config.add_options()
			("bw-epsilon,P", po::value< double >(&Settings::bwEpsilon)->default_value(1e-5),
				"The value within which log-likelihood scores must converge for Baum-Welch.")
			("bw-max-iterations,M", po::value< unsigned int >(&Settings::bwMaxIter)->default_value(0),
				"The maximum number of iterations for Baum-Welch.  A value of zero (0) means no maximum and Baum-Welch is performed to convergence.")
			("bw-pc-during,R", po::value< bool >(&Settings::bwPCDuring)->default_value(true),
				"Specify whether to add pseudocounts during each Baum-Welch iteration (true) or after it is finished (false).")
			("clusters,c", po::value< string >(&clustersFile)->default_value("clusters.txt"),
				"The file from which to read the clusters.  The clusters can have alignment data; this is ignored unless using --prealigned=1.")
			("cpmatrix,x", po::value< string >(&cpFile)->default_value("cp.txt"),
				"The file from which to read the conditional probabilities for a substitution-matrix-based pseudocount method.  A 26x26 matrix is read in, where matrix[a][b] = P(a|b).  In order to use this, we must have --pc-method=2.")
			("favour-m,F", po::value< bool >(&Settings::initFavourM)->default_value(true),
				"Specify whether to initailise the profile HMM with a favouring of the transitions to match states.  Does nothing if --prealigned=1.")
			("frequencies,q", po::value< string >(&qFile)->default_value("q.txt"),
				"The file containing the character frequencies.")
			("freq-scale,Q", po::value< double>(&Settings::initQScale)->default_value(3500),
				"The scaling factor for the frequencies which are used as parameters for a Dirichlet distribution when initializing a random model.")
			("inserts,i", po::value<bool>(&Settings::initInserts)->default_value(true),
				"Specify whether to construct the profile HMM from the clusters using both match and insert states (true) or using match states only (false) when building using --prealigned=1.")
			("pc-method,m", po::value< unsigned int >(&pcMethod)->default_value(2),
				"Specify the pseudocount method used.  Options are 0 for a constant being added to each value (set the weight to the alphabet size for Laplace's rule), 1 for a background-frequency-based pseudocount, and 2 for a substitution-matrix-based pseudocount.  Option 2 requires the --cpmatrix parameter to be set.")
			("pc-weight,w", po::value< double >(&Settings::pcWeight)->default_value(5.0),
				"The weight used for the emission-state pseudocounts.")
			("prealigned,p", po::value< bool>(&prealigned)->default_value(false),
				"Specify whether to initialize the profile HMM using pre-existing alignment information, and then perform Baum-Welch to optimize the model and generate new alignments.")
			("seed,S", po::value< unsigned long >(&Settings::rngSeed)->default_value(time(NULL)),
				"The seed used for the random number generator.  Useful only if --prealigned=0.  If no value is specified, a call to std::time(NULL) is used.")
			("singles,s", po::value<bool>(&singles)->default_value(false),
				"Specify whether to score against single-word clusters.")
			("transition-a,a", po::value< double >(&Settings::initAlpha)->default_value(5),
				"The alpha-values used for a uniform Dirichlet distribution from which the initial transition probabilities are sampled.")
			;

		po::options_description hidden("Hidden options");
		hidden.add_options()
			("word", po::value<string>(&word)->default_value(""), "Specifies the word that is scored against the clusters.")
			;

		po::options_description cmdline_options;
		cmdline_options.add(generic).add(config).add(hidden);

		po::options_description config_file_options;
		config_file_options.add(config).add(hidden);

		po::options_description visible("Allowed options");
		visible.add(generic).add(config);
	    
		po::positional_options_description pos;
		pos.add("word", 1);

		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(pos).run(), vm);

		ifstream ifs("findcluster.cfg");
		store(parse_config_file(ifs, config_file_options), vm);
		po::notify(vm);

		if (vm.count("help"))
		{
			cout << visible << endl << endl;
			cout << "Any and all parameters can be specified via a findcluster.cfg file.  Use the following format:" << endl << endl;
			cout << "#Anything after a # on a line is not processed (i.e. is a comment)" << endl;
			cout << "param=value #e.g. prealigned=1" << endl;

			return 0;
		}

		if (pcMethod == CONSTANT)
			Settings::pcMethod = CONSTANT;
		else if (pcMethod == BACKGROUND)
			Settings::pcMethod = BACKGROUND;
		else if (pcMethod == SUBSTITUTION)
			Settings::pcMethod = SUBSTITUTION;
		else
		{
			cout << "You provided an invalid pseudocount method.  Valid values are:" << endl;
			cout << "0 for a constant value added to each emission count (Laplace's rule is implemented by setting the weight to the alphabet size)" << endl;
			cout << "1 for a background-frequency-based pseudocount" << endl;
			cout << "2 for a substitution-matrix-based pseudocount, which requires the --cpmatrix parameter to be set" << endl;
			
			return 0;
		}

		if (!word.empty())
		{
			for (unsigned int i = 0; i < word.size(); ++i)
				wordInted.push_back(ctoi(word[i]));

			if (readClusters(clustersFile, labels, clusters))
			{
				double q[ALPHABET];

				if (readCharFreqs(qFile, q))
				{
					Settings::q = q;
					vector< string > labelsScored;
					vector< double > scores;

					if (!prealigned)
						cout << "RNG seed: " << Settings::rngSeed << endl << endl;

					if (Settings::pcMethod == SUBSTITUTION)
					{
						double** cp = new double*[ALPHABET];

						for (int i = 0; i < ALPHABET; ++i)
							cp[i] = new double[ALPHABET];

						if (readCPs(cpFile, cp))
							Settings::cp = cp;
						else
						{
							cout << "Couldn't open conditional probability matrix file!" << endl;

							for (int i = 0; i < ALPHABET; ++i)
								delete cp[i];

							delete cp;

							return 1;
						}
					}

					for (unsigned int i = 0; i < clusters.size(); ++i)
						if ((!singles && clusters[i].size() > 1) || singles)
						{
							ProfileHMM *phmm;

							vector< vector< unsigned int > > cluster;
							double avgLength = 0.0;

							for (unsigned int j = 0; j < clusters[i].size(); ++j)
							{
								vector< unsigned int > word;

								for (unsigned int k = 0; k < clusters[i][j].size(); ++k)
									if (clusters[i][j][k] < ALPHABET)
										word.push_back(clusters[i][j][k]);

								cluster.push_back(word);

								avgLength += word.size();
							}

							avgLength /= cluster.size();

							if (prealigned)
							{
								cout << "Building initial profile HMM from pre-aligned cluster \"" << labels[i] << "\"..." << endl;

								phmm = new ProfileHMM(clusters[i]);
							}
							else
							{
								cout << "Generating initial profile HMM..." << endl;

								phmm = new ProfileHMM(dtoi(avgLength));
							}

							cout << "Initial profile HMM built.  Optimizing to cluster \"" << labels[i] << "\"..." << endl;

							cout << "Profile HMM optimized in " << phmm->baumWelch(cluster) << " Baum-Welch iterations." << endl << endl;

							labelsScored.push_back(labels[i]);
							scores.push_back(phmm->logodds(wordInted));

							delete phmm;
						}

					vector< double > scoresSorted(scores);
					sort(scoresSorted.begin(), scoresSorted.end());
					reverse(scoresSorted.begin(), scoresSorted.end());

					printScores(labelsScored, scores);

					bool stop = false;
					for (unsigned int i = 0; i < scores.size() && !stop; ++i)
						if (scores[i] == scoresSorted[0])
						{
							cout << "Highest score was obtained on cluster \"" << labelsScored[i] << "\"." << endl;
							stop = true;
						}

					if (Settings::cp != NULL)
					{
						for (int i = 0; i < ALPHABET; ++i)
							delete Settings::cp[i];

						delete Settings::cp;
					}
				}
				else
				{
					cerr << "Couldn't open frequency data file!" << endl;

					retVal = 1;
				}
			}
			else
			{
				cerr << "Couldn't open clusters file!" << endl;

				retVal = 1;
			}
		}
		else
		{
			cerr << "No word specified!" << endl;

			retVal = 1;
		}
	}
	catch (exception &e)
	{
		cout << "Error: " << e.what() << endl;

		retVal = 1;
	}

	return retVal;
}


void printScores(vector< string > labels, vector< double > scores)
{
	assert(labels.size() == scores.size()); // just in case

	cout << "SCORES" << endl << "------" << endl;

	for (unsigned int i = 0; i < labels.size(); ++i)
		cout << " -" << labels[i] << " " << scores[i] << endl;

	cout << endl;
}
