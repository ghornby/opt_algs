/***

This file is part of the OptAlgs C++ library originally developed
by Gregory Hornby.

This code is released under the BSD-3 license:
  http://opensource.org/licenses/BSD-3-Clause
See the file "license.txt" in the root directory for full details.

***/


/*

 author: Gregory S. Hornby
 file: opt_alg.cpp

 Description:
 The root class for optimization algorithms and defines the API.

*/

#include "opt_alg.h"

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


ostream& operator<<(ostream& os, const OptAlg& a)
{
  return a.write(os);
}

istream& operator>>(istream& is, OptAlg& a)
{
  return a.read(is);
}


const int OPTALG_VERSION = 1;

const int MAX_TRIES = 1;
const int ALG_VERSION = 1;


bool OptAlg :: is_print_log_ = true;
const string OptAlg :: class_name_("OptAlg");


/////////////////////////////////////////////////////////////////////////

vector<Individual*> OptAlg :: individual_types_;
void OptAlg :: add_individual_type(const Individual* ind)
{
  for (vector<Individual*>::size_type i=0; i<individual_types_.size(); i++) {
    if (individual_types_[i]->get_class_name() == ind->get_class_name()) {
      // Already have this class, so quite:
      return;
    }
    /*
    if (OptAlg::is_print_log_) {
      Logger::Msg(2) << "OptAlg :: add_individual_type() - " << ind->get_class_name() << endl;
    }
    */
  }
  individual_types_.push_back(ind->new_instance());
}
const Individual* OptAlg :: get_individual_type(const string& name)
{
  for (vector<Individual*>::size_type i=0; i<individual_types_.size(); i++) {
    if (individual_types_[i]->get_class_name() == name) {
      return individual_types_[i];
    }
  }
  return 0;
}
void OptAlg :: clear_individual_types()
{
  for (vector<Individual*>::size_type i=0; i<individual_types_.size(); i++) {
    delete individual_types_[i];
    individual_types_[i] = 0;
  }
  individual_types_.clear();
}


///////

vector<OptAlg*> OptAlg :: alg_types_;
void OptAlg :: add_alg_type(const OptAlg* alg)
{
  for (vector<OptAlg*>::size_type i=0; i<alg_types_.size(); i++) {
    if (alg_types_[i]->get_class_name() == alg->get_class_name()) {
      // Already have this class, so quite:
      return;
    }
  }
  alg_types_.push_back(alg->new_instance());
}
const OptAlg* OptAlg :: get_alg_type(const string& name)
{
  for (vector<OptAlg*>::size_type i=0; i<alg_types_.size(); i++) {
    if (alg_types_[i]->get_class_name() == name) {
      return alg_types_[i];
    }
  }
  return 0;
}
void OptAlg :: clear_alg_types()
{
  for (vector<OptAlg*>::size_type i=0; i<alg_types_.size(); i++) {
    delete alg_types_[i];
    alg_types_[i] = 0;
  }
  alg_types_.clear();
}


void OptAlg :: clear_types()
{
  OptAlg::clear_individual_types();
  OptAlg::clear_alg_types();
}


void OptAlg :: set_random_seed(unsigned long seed)
{
  CmnClass::seed_random(seed);
}



/********************************************************/
/********************************************************/
/********************************************************/

/***** Initialization functions *****/


/**
 * Main constructor for optimization algorithm
 *
 *
 */
OptAlg :: OptAlg(const FitnessFunc* func, const Individual* ind_template)
  : is_maximizing_(false), num_evaluations_(0), ind_template_(0),
    ind_seed_(0)
{
  is_print_debug_info_ = false;
  is_run_sample_ = false;
  is_save_log_ = false;
  is_save_log_ = false;


  if (ind_template) {
    ind_template_ = ind_template->new_instance();
    ind_template_->duplicate_settings(ind_template);
  } else {
    ind_template_ = 0;
  }

  fitness_func_ = 0;
  if (func) {
    fitness_func_ = func->new_instance();
  }
}



OptAlg :: ~OptAlg()
{
  if (ind_template_) {
    delete ind_template_;
    ind_template_ = 0;
  }

  if (ind_seed_) {
    delete ind_seed_;
    ind_seed_ = 0;
  }

  if (fitness_func_) {
    delete fitness_func_;
    fitness_func_ = 0;
  }
}



void OptAlg :: copy_settings(const OptAlg* p_src)
{
  is_maximizing_ = p_src->is_maximizing_;

  is_run_sample_ = p_src->is_run_sample_;
  is_save_best_ = p_src->is_save_best_;
  is_save_log_ = p_src->is_save_log_;
  log_file_name_ = p_src->log_file_name_;

  is_print_debug_info_ = p_src->is_print_debug_info_;

  num_evaluations_ = p_src->num_evaluations_;

  ind_template_ = p_src->ind_template_->new_instance();
  ind_template_->duplicate(p_src->ind_template_);

  ind_seed_ = 0;
  if (p_src->ind_seed_) {
    ind_seed_ = p_src->ind_template_->new_instance();
    ind_seed_->duplicate(p_src->ind_seed_);
  }
}



void OptAlg :: clear()
{
  is_print_debug_info_ = false;
  is_run_sample_ = false;
  is_save_log_ = false;
  is_save_log_ = false;

  if (ind_template_) {
    ind_template_->clear();
  }
  if (ind_seed_) {
    ind_seed_->clear();
    delete ind_seed_;
    ind_seed_ = 0;
  }

  num_evaluations_ = 0;
}

void OptAlg :: init_variables()
{
  num_evaluations_ = 0;
}


int OptAlg :: set_logfile(const char *n)
{
  if (n == 0) {
    return false;
  }

  if (n[0] == 0) {
    return false;
  }

  log_file_name_.assign(n);
  log_file_name_ += ".log";

  is_save_log_ = true;

  // Open log file.
  log_file_.open(log_file_name_.c_str());
  if (is_print_log_) {
    Logger::Msg(2) << "Opening log: " << log_file_name_ << "\n";
  }

  return true;
}


  //int OptAlg :: configure(const char *fname, int verbose)
int OptAlg :: configure(const char *, int)
{
  Logger::ErrorDTS() << "opt_alg :: configure() - no longer supported.  while (1) ;" << endl;
  CmnClass::Logger::flush_log_error();
  while (1) ;

  return true;
}




// ************************************************************ //


void OptAlg :: set_print_debug(bool f)
{
  is_print_debug_info_ = f;
}

void OptAlg :: set_save_best(bool b)
{
  is_save_best_ = b;
}



/*********************************************************************************/
/*********************************************************************************/


unsigned int OptAlg :: get_num_evals() const
{
  return num_evaluations_;
}

void OptAlg :: increment_num_evals(unsigned int num)
{
  num_evaluations_ += num;

  if (((num_evaluations_ < 101) && (num_evaluations_ % 10) == 0) ||
      ((num_evaluations_ < 1001) && (num_evaluations_ % 100) == 0) ||
      ((num_evaluations_ < 10010) && (num_evaluations_ % 1000) == 0) ||
      ((num_evaluations_ < 100100) && (num_evaluations_ % 5000) == 0) ||
      ((num_evaluations_ < 1001000) && (num_evaluations_ % 50000) == 0) ||
      ((num_evaluations_ % 100000) == 0)) {
    cout << num_evaluations_ << " " << fitness_func_->get_best_fitness() << endl;
  }
}
void OptAlg :: increment_num_evals()
{
  increment_num_evals(1);
}


/*********************************************************************************/
/*********************************************************************************/


// Note: probably need to give this function an existing
//  individual to read into.
// For more ellaborate individual types, this individual
//  would already be properly configured with a translator function.
bool OptAlg :: read_sample_individ(const char *fname)
{
  if (is_run_sample_) {
    printf("opt_alg :: read_sample_individ() - already loaded individual.\n");
    return false;
  }

  if (ind_template_ == 0) {
 //   ind_template_ = new Individual();
  }

  return ind_template_->read(fname);
}


bool OptAlg :: read_seed_individ(const char *fname)
{
  if (ind_template_ == 0) {
 //   ind_seed_ = new Individual();
  }

  return ind_seed_->read(fname);
}


// ************************************************************** //
// ************************************************************** //


istream& OptAlg :: read(istream& istr)
{
  string line;
  vector<string> words;

  getline(istr, line);
  CmnClass::split(line, ' ', words);
  if (words.size() != 2) {
    throw runtime_error("opt_alg :: read() - error with initial line 1.");
  }

  if (words[0] != "OptAlgV") {
    throw runtime_error("opt_alg :: read() - error, expecting OptAlgV, got:" + line);
  }

  int version = CmnClass::from_string<int>(words[1]);
  if (version != OPTALG_VERSION) {
    throw runtime_error("opt_alg :: read() - error, bad version:" + words[1]);
  }  


  getline(istr, line);
  CmnClass::split(line, ' ', words);
  if (words.size() != 1) {
    throw runtime_error("opt_alg :: read() - error reading log_file_name:" + line);
  }
  log_file_name_ = words[0];
  if (log_file_name_ == "-") {
    log_file_name_.clear();
  }



  getline(istr, line);
  CmnClass::split(line, ' ', words);
  if (words.size() != 5) {
    Logger::Msg(2) << "opt_alg :: read() - error parsing: ]" << line << "[" << endl;
    throw runtime_error("opt_alg :: read() - error config line: " + line);
  }
  is_maximizing_ = CmnClass::from_string<bool>(words[1]);
  is_run_sample_ = CmnClass::from_string<bool>(words[2]);
  is_save_best_ = CmnClass::from_string<bool>(words[3]);
  is_save_log_ = CmnClass::from_string<bool>(words[4]);


  getline(istr, line);
  CmnClass::split(line, ' ', words);
  if (words.size() != 1) {
    Logger::Msg(2) << "opt_alg :: read() - error parsing: ]" << line << "[" << endl;
    throw runtime_error("opt_alg :: read() - error config line: " + line);
  }
  num_evaluations_ = CmnClass::from_string<int>(words[0]);

  return istr;
}


bool OptAlg :: read(const char *fname)
{
  ifstream file(fname);
  try {
    read(file);
    //  } catch (OptAlgError e) {
  } catch (...) {
    return false;
  }

  file.close();

  return true;
}


ostream& OptAlg :: write_header(ostream& ostr) const
{
  ostr << "# OptAlgV " << OPTALG_VERSION << endl;

  ostr << "# ";
  if (log_file_name_[0] != 0) {
    ostr << log_file_name_ << endl;
  } else {
    ostr << "-" << endl;
  }

  ostr << "# ";
  ostr << is_maximizing_ << " "    
       << is_run_sample_ << " "
       << is_save_best_ << " "
       << is_save_log_ << endl;

  ostr << "# ";
  ostr << num_evaluations_ << endl;

  return ostr;
}


ostream& OptAlg :: write_template(ostream& ostr) const
{
  ostr << "OptAlgV " << OPTALG_VERSION << endl;

  return ostr;
}


ostream& OptAlg :: write(ostream& ostr) const
{
  ostr << "OptAlgV " << OPTALG_VERSION << endl;

  if (log_file_name_.size() != 0) {
    ostr << log_file_name_ << endl;
  } else {
    ostr << "-" << endl;
  }

  ostr << is_maximizing_ << " "    
       << is_run_sample_ << " "
       << is_save_best_ << " "
       << is_save_log_ << endl;

  ostr << num_evaluations_ << endl;

  return ostr;
}


bool OptAlg :: write(const char *fname) const
{
  ofstream file(fname);
  write(fname);
  file.close();

  return true;
}


void OptAlg :: make_filename(const char *ending, char *fname) const
{
  int i;

  for (i=0; i<40; i++) {
    fname[i] = class_name_[i];
    if (fname[i] == 0)
      break;
  }

  strcpy(&fname[i], ending);
  /*
  if (num_runs_ > 1) {
    sprintf(&fname[i], ".%d%s", run_number_, ending);
  } else {
    strcpy(&fname[i], ending);
  }
  */
}


bool OptAlg :: write_backup(void) const
{
  char fname[60];

  make_filename("_bk.pop", fname);
  return write(fname);
}

