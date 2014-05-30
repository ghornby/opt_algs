Project Title: Opt Algs
Auther: Gregory Hornby
====================================

Description:

This repository contains C++ code for a library of optimization
algorithms.  My goal is to implement some of the ones described
in Spall's book, "Introduction to Stochastic Search and Optimization"
as well as add an evolutionary algorithm or two.  Also, I would
then like to hybridize these with my Age-Layered Population
Structure (ALPS) meta-heuristic.


Directory structure is as follows:
o examples: this contains example source code for using this library.
o opt_algs.X.Y: contains the source code for version X.Y.


At some point I'm likely to put these functions in a namespace,
such as OptAlgs.


The Makefile has options for installing the library in your
user directory using, 'make installusr' (which assumes that you
already have the directories: ~/include and ~/lib) or to install
it globally on your computer with, 'make installglobal'.

These installs will put copies of files in the following directories:

installusr:
  headers: ~/include/optalg.X.Y
  library: ~/lib/optalg.X.Y.a

installglobal:
  headers: /usr/local/include/optalg.X.Y
  library: /usr/local/lib/optalg.X.Y.a


In addition, links will be created to these as simply:
   optalg.X or optalg.X.a

New versions of the common class library which keep the
same API will increment the Y counter.  That is,
  optalg.1.0 => optalg.1.1

When/if a rewrite occurs in which the API changes, then X
will be incremented:
  optalg.1.0 => optalg.2.0



Dependencies:

This library depends on code from my "common_classes" library.



==========================================================

License:

This project is being developed in my free time and is released
under the BSD-3 license:
  http://opensource.org/licenses/BSD-3-Clause

See the file ./license.txt for details.



