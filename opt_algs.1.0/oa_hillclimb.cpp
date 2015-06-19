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
 The source code for a simple hill-climbing search algorithm.

*/

#include "oa_hillclimb.h"

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


ostream& operator<<(ostream& os, const HillClimber& a)
{
  return a.write(os);
}

istream& operator>>(istream& is, HillClimber& a)
{
  return a.read(is);
}


const int OA_RANDOM_VERSION = 1;

const int MAX_TRIES = 1;
const int ALG_VERSION = 1;


const string HillClimber :: class_name_("HillClimber");



/*********************************************************************************/
/*********************************************************************************/

/***** Initialization functions *****/


/**
 * Main constructor for optimization algorithm
 *
 *
 */
HillClimber :: HillClimber(const FitnessFunc* func, const Individual* ind_template)
  : OptAlg(func, ind_template), individ_best_(0), individ_new_(0)
{
  if (ind_template) {
    individ_best_ = ind_template->new_instance();
    individ_best_->duplicate(ind_template);

    individ_new_ = ind_template->new_instance();
  }
  /*
  CmnClass::Logger::Msg() << "individ_best #evals: "
			  << individ_best_->get_num_evaluations() << endl;
  CmnClass::Logger::Msg() << "individ_template #evals: "
			  << ind_template->get_num_evaluations() << endl;
  */
}



HillClimber :: ~HillClimber()
{
}



void HillClimber :: copy_settings(const HillClimber* p_src)
{
  OptAlg::copy_settings(p_src);
}


OptAlg* HillClimber :: new_instance() const
{
  return new HillClimber(fitness_func_, ind_template_);
}

OptAlg* HillClimber :: make_copy() const
{
  if (is_print_log_) {
    Logger::Msg(2) << "oa_random :: make_copy()" << endl;
  }
  
  // Returns a copy of this.
  OptAlg* c = new_instance();
  c->copy_settings(this);

  return c;
}


void HillClimber :: clear()
{
  OptAlg::clear();
}

void HillClimber :: init_variables()
{
  OptAlg::init_variables();
}


  //int HillClimber :: configure(const char *fname, int verbose)
int HillClimber :: configure(const char *, int)
{
  Logger::ErrorDTS() << "oa_random :: configure() - no longer supported.  while (1) ;" << endl;
  CmnClass::Logger::flush_log_error();
  while (1) ;

  return true;
}



/*********************************************************************************/
/*********************************************************************************/

/*   Search/Optimization functions     */



void HillClimber :: do_search_step()
{
  //  individ_best_->write("best.ind");

  if (individ_best_->is_keep(is_maximizing_)) {
    // If best individual is good enough to keep,
    // copy and mutate:
    //    CmnClass::Logger::Msg() << "HC :: do_search_step() - mutating" << endl;
    individ_new_->duplicate(individ_best_);
    //    individ_best_->write("duplicate.ind");
    individ_new_->mutate();

  } else {
    //    CmnClass::Logger::Msg() << "HC :: do_search_step() - make random" << endl;
    // Best individual is flawed (typically at start of run)
    // so create a new random individual to start from:
    individ_new_->make_random();
  }

  //  individ_best_->write("new.ind");

  double val = fitness_func_->evaluate(individ_new_);

  if (individ_best_->compare_fitness(is_maximizing_, individ_new_) == -1) {
    // individ_new_ is better, so make individ_best_ a copy of it.
    individ_best_->duplicate(individ_new_);
  }

  increment_num_evals();

  /*
  unsigned int num_evals = get_num_evals();
  stringstream fname;
  fname << "step_";
  if (num_evals < 10) {
    fname << "0";
  }
  if (num_evals < 100) {
    fname << "0";
  }
  if (num_evals < 1000) {
    fname << "0";
  }
  fname << num_evals << ".ind";

  string fn_str(fname.str());
  individ_new_->write(fn_str.c_str());
  */
}



/*********************************************************************************/
/*********************************************************************************/

/*    I/O functions    */


void HillClimber :: print_current_individs(std::ostream& ostr)
{
  if (individ_best_) {
    ostr << "Individ best:" << endl;
    ostr << individ_best_ << endl;
  } else {
    ostr << "No current individs." << endl;
  }
}
void HillClimber :: print_current_individs()
{
  print_current_individs(cout);
}



istream& HillClimber :: read(istream& istr)
{
  string line;
  vector<string> words;

  // Read version:
  getline(istr, line);
  CmnClass::split(line, ' ', words);
  if (words.size() != 2) {
    throw runtime_error("oa_random :: read() - error with initial line 1.");
  }

  if (words[0] != "HillClimberV") {
    throw runtime_error("oa_random :: read() - error, expecting HillClimberV, got:" + line);
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


bool HillClimber :: read(const char *fname)
{
  ifstream file(fname);
  try {
    read(file);
    //  } catch (HillClimberError e) {
  } catch (...) {
    return false;
  }

  file.close();

  return true;
}


ostream& HillClimber :: write_header(ostream& ostr) const
{
  ostr << "# HillClimberV " << OA_RANDOM_VERSION << endl;

  OptAlg::write_header(ostr);


  return ostr;
}


ostream& HillClimber :: write_template(ostream& ostr) const
{
  ostr << "HillClimberV " << OA_RANDOM_VERSION << endl;

  return ostr;
}


ostream& HillClimber :: write(ostream& ostr) const
{
  ostr << "HillClimberV " << OA_RANDOM_VERSION << endl;

  OptAlg::write(ostr);

  return ostr;
}


bool HillClimber :: write(const char *fname) const
{
  ofstream file(fname);
  write(fname);
  file.close();

  return true;
}



bool HillClimber :: write_backup(void) const
{
  char fname[60];

  make_filename("_bk.pop", fname);
  return write(fname);
}

