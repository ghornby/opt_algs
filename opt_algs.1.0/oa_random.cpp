/***

This file is part of the OptAlgs C++ library originally developed
by Gregory Hornby.

This code is released under the BSD-3 license:
  http://opensource.org/licenses/BSD-3-Clause
See the file "license.txt" in the root directory for full details.

***/


/*

 author: Gregory S. Hornby
 file: oa_random.cpp

 Description:
 The source code for a simple random-search class.

*/

#include "oa_random.h"

// Project-specific includes:
#include "individual.h"
#include "oa_function.h"


// Developed includes, generic:
#include "cmn_logger.h"
#include "cmn_random.h"
#include "cmn_string_proc.h"

#include "cmn_logger.h"
using CmnClass::Logger;



// Standard system includes:
#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace std;


ostream& operator<<(ostream& os, const OA_Random& a)
{
  return a.write(os);
}

istream& operator>>(istream& is, OA_Random& a)
{
  return a.read(is);
}


const int OA_RANDOM_VERSION = 1;

const int MAX_TRIES = 1;
const int ALG_VERSION = 1;


const string OA_Random :: class_name_("OA_Random");



/*********************************************************************************/
/*********************************************************************************/

/***** Initialization functions *****/


/**
 * Main constructor for optimization algorithm
 *
 *
 */
OA_Random :: OA_Random(const FitnessFunc* func, const Individual* ind_template)
  : OptAlg(func, ind_template), individ_best_(0), individ_new_(0)
{
  individ_best_ = ind_template->new_instance();
  individ_new_ = ind_template->new_instance();
}



OA_Random :: ~OA_Random()
{
}



void OA_Random :: copy_settings(const OA_Random* p_src)
{
  OptAlg::copy_settings(p_src);
}


OptAlg* OA_Random :: new_instance() const
{
  return new OA_Random(fitness_func_, ind_template_);
}

OptAlg* OA_Random :: make_copy() const
{
  if (is_print_log_) {
    Logger::Msg(2) << "oa_random :: make_copy()" << endl;
  }
  
  // Returns a copy of this.
  OptAlg* c = new_instance();
  c->copy_settings(this);

  return c;
}


void OA_Random :: clear()
{
  OptAlg::clear();
}

void OA_Random :: init_variables()
{
  OptAlg::init_variables();
}


  //int OA_Random :: configure(const char *fname, int verbose)
int OA_Random :: configure(const char *, int)
{
  Logger::ErrorDTS() << "oa_random :: configure() - no longer supported.  while (1) ;" << endl;
  CmnClass::Logger::flush_log_error();
  while (1) ;

  return true;
}



/*********************************************************************************/
/*********************************************************************************/

/*   Search/Optimization functions     */



void OA_Random :: do_search_step()
{
  individ_new_->make_random();
  double val = fitness_func_->evaluate(individ_new_);

  if (individ_best_->compare_fitness(is_maximizing_, individ_new_) == -1) {
    // individ_new_ is better, so make individ_best_ a copy of it.
    individ_best_->duplicate(individ_new_);
  }

  increment_num_evals();
}



/*********************************************************************************/
/*********************************************************************************/

/*    I/O functions    */


void OA_Random :: print_current_individs(std::ostream& ostr)
{
  if (individ_best_) {
    ostr << "Individ best:" << endl;
    ostr << individ_best_ << endl;
  } else {
    ostr << "No current individs." << endl;
  }
}
void OA_Random :: print_current_individs()
{
  print_current_individs(cout);
}



istream& OA_Random :: read(istream& istr)
{
  string line;
  vector<string> words;

  // Read version:
  getline(istr, line);
  CmnClass::split(line, ' ', words);
  if (words.size() != 2) {
    throw runtime_error("oa_random :: read() - error with initial line 1.");
  }

  if (words[0] != "OA_RandomV") {
    throw runtime_error("oa_random :: read() - error, expecting OA_RandomV, got:" + line);
  }

  int version = CmnClass::from_string<int>(words[1]);
  if (version != OA_RANDOM_VERSION) {
    throw runtime_error("oa_random :: read() - error, bad version:" + words[1]);
  }  


  // Have parent class read whatever it stores:
  OptAlg::read(istr);


  // Read state/things for this class:


  return istr;
}


bool OA_Random :: read(const char *fname)
{
  ifstream file(fname);
  try {
    read(file);
    //  } catch (OA_RandomError e) {
  } catch (...) {
    return false;
  }

  file.close();

  return true;
}


ostream& OA_Random :: write_header(ostream& ostr) const
{
  ostr << "# OA_RandomV " << OA_RANDOM_VERSION << endl;

  OptAlg::write_header(ostr);


  return ostr;
}


ostream& OA_Random :: write_template(ostream& ostr) const
{
  ostr << "OA_RandomV " << OA_RANDOM_VERSION << endl;

  return ostr;
}


ostream& OA_Random :: write(ostream& ostr) const
{
  ostr << "OA_RandomV " << OA_RANDOM_VERSION << endl;

  OptAlg::write(ostr);

  return ostr;
}


bool OA_Random :: write(const char *fname) const
{
  ofstream file(fname);
  write(fname);
  file.close();

  return true;
}



bool OA_Random :: write_backup(void) const
{
  char fname[60];

  make_filename("_bk.pop", fname);
  return write(fname);
}

