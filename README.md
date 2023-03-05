# Profile HMM tools

This repository contains an implementation of Profile HMMs (hidden Markov models) intended for aligning multiple words together as described in our paper, [Multiple Word Alignment with Profile Hidden Markov Models](https://aclanthology.org/N09-3008/). Profile HMMs are used in biological sequence analysis, specifically for the alignment of sequences and the testing of sequences for membership in a family. At the time of implementation, the existing Profile HMM libraries ([SAM](http://web.archive.org/web/20220121114141/https://compbio.soe.ucsc.edu/sam.html) and [HMMER](http://hmmer.org/)) were too specific to biological sequences and difficult to adapt to word-related tasks.

## Using the tools or the library in your work

The tools and the library are licensed under the GPL for non-commercial use only. If you use it in your research, please cite the corresponding paper:
* [Aditya Bhargava](https://www.cs.toronto.edu/~aditya/) and [Grzegorz Kondrak](https://webdocs.cs.ualberta.ca/~kondrak/). 2009. [Multiple Word Alignment with Profile Hidden Markov Models](https://aclanthology.org/N09-3008/).  In *Proceedings of Human Language Technologies: The 2009 Annual Conference of the North American Chapter of the Association for Computational Linguistics, Companion Volume: Student Research Workshop and Doctoral Consortium*, pages 43–48, Boulder, Colorado.  Association for Computational Linguistics.

## Building

### Prerequisites

Only Linux is supported.[^windows] Other than the usual build dependencies (like `g++`, `make`, etc.), you will need:
* The [GNU Scientific Library](https://www.gnu.org/software/gsl/)
* The [Boost libraries](https://www.boost.org/)
Most (if not all) distributions provide packages for these; you'll have to consult your distribution's package repository to determine what they're called. When installing the prerequisites, make sure you include the development header packages as well as the runtime libraries (most, but I think not all, distributions separate the two).

### Compilation

Given that the code was written around 2008, it's somewhat miraculous that it still compiles with just a few warnings (that I am choosing to ignore). Simply check out the code from this repository and then run `make` from within the `phmmtools` subdirectory.[^compiles]

## Usage

The code is split into a Profile HMM library and a pair of CLI (command-line interface) tools.

### Profile HMM library

The library, `libprofilehmm`, provides all of the necessary functions to create and train a Profile HMM, as well as align or score words to a model. The library is undocumented, though, so the best way to make use of the library by itself is to take a look at the source code for the CLI tools.

### Profile HMM tools

We tested Profile HMMs for applicability to multiple word alignment and cognate set matching (I refer to "cognate sets" as "clusters" here). There are thus two programs:
1. <code>multialign</code> for multiple alignment; and 
2. <code>findcluster</code> to match a word to a cluster.

For either tool, call the program with the `-h` argument to get a list of options:

#### `multialign`

```
$ ./multialign -h
Allowed options:

Generic options:
  -h [ --help ]                         Prints help message.

Configuration options:
  -o [ --alignments ] arg (=alignments.txt)
                                        The file to which to write the 
                                        alignments.  Note that the file will be
                                        overwritten if it exists.
  -P [ --bw-epsilon ] arg (=1.0000000000000001e-05)
                                        The value within which log-likelihood 
                                        scores must converge for Baum-Welch.
  -M [ --bw-max-iterations ] arg (=0)   The maximum number of iterations for 
                                        Baum-Welch.  A value of zero (0) means 
                                        no maximum and Baum-Welch is performed 
                                        to convergence.
  -R [ --bw-pc-during ] arg (=1)        Specify whether to add pseudocounts 
                                        during each Baum-Welch iteration (true)
                                        or after it is finished (false).
  -c [ --clusters ] arg (=clusters.txt) The file from which to read the 
                                        clusters.  The clusters can have 
                                        alignment data; this is ignored unless 
                                        using --prealigned=1.
  -x [ --cpmatrix ] arg (=cp.txt)       The file from which to read the 
                                        conditional probabilities for a 
                                        substitution-matrix-based pseudocount 
                                        method.  A 26x26 matrix is read in, 
                                        where matrix[a][b] = P(a|b).  In order 
                                        to use this, we must have 
                                        --pc-method=2.
  -F [ --favour-m ] arg (=1)            Specify whether to initailise the 
                                        profile HMM with a favouring of the 
                                        transitions to match states.  Does 
                                        nothing if --prealigned=1.
  -q [ --frequencies ] arg (=q.txt)     The file containing the character 
                                        frequencies.
  -Q [ --freq-scale ] arg (=3500)       The scaling factor for the frequencies 
                                        which are used as parameters for a 
                                        Dirichlet distribution when 
                                        initializing a random model.
  -i [ --inserts ] arg (=1)             Specify whether to construct the 
                                        profile HMM from the clusters using 
                                        both match and insert states (true) or 
                                        using match states only (false) when 
                                        building using --prealigned=1.
  -m [ --pc-method ] arg (=2)           Specify the pseudocount method used.  
                                        Options are 0 for a constant being 
                                        added to each value (set the weight to 
                                        the alphabet size for Laplace's rule), 
                                        1 for a background-frequency-based 
                                        pseudocount, and 2 for a 
                                        substitution-matrix-based pseudocount. 
                                        Option 2 requires the --cpmatrix 
                                        parameter to be set.
  -w [ --pc-weight ] arg (=5)           The weight used for the emission-state 
                                        pseudocounts.
  -p [ --prealigned ] arg (=0)          Specify whether to initialize the 
                                        profile HMM using pre-existing 
                                        alignment information, and then perform
                                        Baum-Welch to optimize the model and 
                                        generate new alignments.
  -S [ --seed ] arg (=1678045268)       The seed used for the random number 
                                        generator.  Useful only if 
                                        --prealigned=0.  If no value is 
                                        specified, a call to std::time(NULL) is
                                        used.
  -a [ --transition-a ] arg (=5)        The alpha-values used for a uniform 
                                        Dirichlet distribution from which the 
                                        initial transition probabilities are 
                                        sampled.


Any and all parameters can be specified via a multialign.cfg file.  Use the following format:

#Anything after a # on a line is not processed (i.e. is a comment)
param=value #e.g. prealigned=1
```

#### `findcluster`

```
$ ./findcluster -h
Allowed options:

Generic options:
  -h [ --help ]                         Prints help message.

Configuration options:
  -P [ --bw-epsilon ] arg (=1.0000000000000001e-05)
                                        The value within which log-likelihood 
                                        scores must converge for Baum-Welch.
  -M [ --bw-max-iterations ] arg (=0)   The maximum number of iterations for 
                                        Baum-Welch.  A value of zero (0) means 
                                        no maximum and Baum-Welch is performed 
                                        to convergence.
  -R [ --bw-pc-during ] arg (=1)        Specify whether to add pseudocounts 
                                        during each Baum-Welch iteration (true)
                                        or after it is finished (false).
  -c [ --clusters ] arg (=clusters.txt) The file from which to read the 
                                        clusters.  The clusters can have 
                                        alignment data; this is ignored unless 
                                        using --prealigned=1.
  -x [ --cpmatrix ] arg (=cp.txt)       The file from which to read the 
                                        conditional probabilities for a 
                                        substitution-matrix-based pseudocount 
                                        method.  A 26x26 matrix is read in, 
                                        where matrix[a][b] = P(a|b).  In order 
                                        to use this, we must have 
                                        --pc-method=2.
  -F [ --favour-m ] arg (=1)            Specify whether to initailise the 
                                        profile HMM with a favouring of the 
                                        transitions to match states.  Does 
                                        nothing if --prealigned=1.
  -q [ --frequencies ] arg (=q.txt)     The file containing the character 
                                        frequencies.
  -Q [ --freq-scale ] arg (=3500)       The scaling factor for the frequencies 
                                        which are used as parameters for a 
                                        Dirichlet distribution when 
                                        initializing a random model.
  -i [ --inserts ] arg (=1)             Specify whether to construct the 
                                        profile HMM from the clusters using 
                                        both match and insert states (true) or 
                                        using match states only (false) when 
                                        building using --prealigned=1.
  -m [ --pc-method ] arg (=2)           Specify the pseudocount method used.  
                                        Options are 0 for a constant being 
                                        added to each value (set the weight to 
                                        the alphabet size for Laplace's rule), 
                                        1 for a background-frequency-based 
                                        pseudocount, and 2 for a 
                                        substitution-matrix-based pseudocount. 
                                        Option 2 requires the --cpmatrix 
                                        parameter to be set.
  -w [ --pc-weight ] arg (=5)           The weight used for the emission-state 
                                        pseudocounts.
  -p [ --prealigned ] arg (=0)          Specify whether to initialize the 
                                        profile HMM using pre-existing 
                                        alignment information, and then perform
                                        Baum-Welch to optimize the model and 
                                        generate new alignments.
  -S [ --seed ] arg (=1678045341)       The seed used for the random number 
                                        generator.  Useful only if 
                                        --prealigned=0.  If no value is 
                                        specified, a call to std::time(NULL) is
                                        used.
  -s [ --singles ] arg (=0)             Specify whether to score against 
                                        single-word clusters.
  -a [ --transition-a ] arg (=5)        The alpha-values used for a uniform 
                                        Dirichlet distribution from which the 
                                        initial transition probabilities are 
                                        sampled.


Any and all parameters can be specified via a findcluster.cfg file.  Use the following format:

#Anything after a # on a line is not processed (i.e. is a comment)
param=value #e.g. prealigned=1
```

#### File formats and sample files

##### Clusters

The data format uses, for each cluster, a single line starting with a hash (`#`)[^hash]; anything after the hash is considered a label for the cluster. The next lines in the file are all considered part of the cluster until the next hash is read. Anything past a double-hash (`##`) is not read at all. See [testclusters.txt](examples/testclusters.txt) for an example, which shows a couple of test clusters that I used.

##### Background frequencies

The background frequency file expects 26 lines, each containing the background frequency (i.e. frequency of occurrence) of the corresponding character (A-Z). See [q.txt](examples/q.txt) for an example, which is the one I used in my experiments.

##### Conditional probabilities

The conditional probability file for the substitution matrix pseudocount method expects a 26×26 matrix. Each line corresponds to a character, and each entry on the line is separated by a space. Entry $b$ on line $a$ represents $P(a|b)$. See [cp.txt](examples/cp.txt) for an example, which is the one I used in my experiments.

## TODOs

None of these is likely to be implemented (by me), but I had the following list of TODO items for a v1.0 release:
* Parameterize the alphabet size (currently hard-coded to 26)
* Maximum a posteriori (MAP) model construction; there is some code present that does not work correctly
* Any other TODOs marked in the code
* Add some proper versioning info
* Code profiling to find possible speed-ups
* Documentation

[^windows]: The code used to work just fine on Windows, as all the prerequisites are cross-platform and the code was written to work with both `g++` and MSVC. But I don't have access to (or the patience for) a Windows machine anymore, so I have no idea the extent to which this holds true 15 years after the code was written.
[^compiles]: The code compiles, but whether or not it *works*... 🤷‍♂️
[^hash]: No, the `#` symbol isn't called a "hashtag."[^hashtag]
[^hashtag]: Yes, I'm that kind of pedant.[^descriptivism]
[^descriptivism]: But I'm also a descriptivist, so, you know, no (overt) judgement or anything.
