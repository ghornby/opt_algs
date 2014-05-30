/***

This file is part of the OptAlgs C++ library originally developed
by Gregory Hornby.

This code is released under the BSD-3 license:
  http://opensource.org/licenses/BSD-3-Clause
See the file "license.txt" in the root directory for full details.

***/


/*******************************************************

  File: ind_history.cpp

********************************************************/


#include "ind_history.h"

#include "cmn_string_proc.h"

#include <iostream>
#include <string>
#include <fstream>
#include <stdexcept>


using namespace std;
using namespace CmnClass;


const int Individ_Init_Age = 1;


ostream& operator<<(ostream& os, const IndHistory& hist)
{
  return hist.write(os);
}


istream& operator>>(istream& is, IndHistory& hist)
{
  return hist.read(is);
}



IndHistory :: IndHistory()
  : age_(Individ_Init_Age), how_created_(CREATE_INIT),
    last_op_(0), num_mutate_(0), num_recombine_(0),
    op_count_(1, 0), complexity_(1, 0.0),
    fitness_(1, 0.0), fitness_prev_(1, 0.0)
{
  /*
  op_count_.resize(1, 0);
  fitness_.resize(1, 0.0);
  fitness_prev_.resize(1, 0.0);
  complexity_.resize(1, 0.0);
  */
}

IndHistory :: ~IndHistory()
{

}


void IndHistory :: clear()
{
  age_ = Individ_Init_Age;
  how_created_ = CREATE_INIT;
  num_mutate_ = 0;
  num_recombine_ = 0;

  /*
  op_count_.resize(0);
  fitness_.resize(0);
  fitness_prev_.resize(0);
  complexity_.resize(0);  
  */

  op_count_.clear();
  op_count_.resize(1, 0);

  fitness_.clear();
  fitness_.resize(1, 0.0);

  fitness_prev_.clear();
  fitness_prev_.resize(1, 0.0);

  complexity_.clear();
  complexity_.resize(1, 0.0);
}



void IndHistory :: backup_fitness()
{
  fitness_prev_ = fitness_;
}

void IndHistory :: avg_fitness(const vector<double>& fvec)
{
  if (fitness_.size() != fvec.size()) {
    throw runtime_error("history :: avg_fitness() - error, fitness vectors not same size.");
  }

  for (vector<double>::size_type i = 0; i<fitness_.size(); i++) {
    fitness_[i] = (fitness_[i] + fvec[i]) / 2.0;
  }
}

void IndHistory :: set_fitness(const vector<double>& fvec)
{
  fitness_ = fvec;
}



double IndHistory :: get_fitness() const
{
  return fitness_[0];
}

double IndHistory :: get_fitness(int index) const
{
  return fitness_[index];
}

const vector<double>& IndHistory :: get_fitness_vec() const
{
  return fitness_;
}


double IndHistory :: get_fitness_prev() const
{
  return fitness_prev_[0];
}


double IndHistory :: get_complexity(int index) const
{
  return complexity_[index];
}


int IndHistory :: get_how_created() const
{
  return how_created_;
}

void IndHistory :: set_how_created(int h)
{
  how_created_ = h;
}


int IndHistory :: get_age() const
{
  return age_;
}

int IndHistory :: get_age(int evals, int denom) const
{
  //  return (evals - age_)/denom + Individ_Init_Age;
  return (age_+evals)/denom + Individ_Init_Age;
}

void IndHistory :: set_age(int a)
{
  age_ = a;
}

void IndHistory :: incr_age()
{
  age_++;
}

void IndHistory :: take_older(const IndHistory& h2)
{
  if (h2.age_ > age_) {
    age_ = h2.age_;
  }
}

void IndHistory :: take_younger(const IndHistory& h2)
{
  if (h2.age_ < age_) {
    age_ = h2.age_;
  }
}


void IndHistory :: take_average(const IndHistory& h2)
{
  age_ = (age_ + h2.age_) / 2;
}


void IndHistory :: incr_num_mutate()
{
  num_mutate_++;
}

void IndHistory :: incr_num_recombine()
{
  num_recombine_++;
}



vector<int>* IndHistory :: get_op_count()
{
  return &op_count_;
}


void IndHistory :: print_created_orig(ostream& ofs)
{
  if (how_created_ == CREATE_RANDOM) {
    ofs << "[*:";
  } else if (how_created_ == CREATE_INIT) {
    ofs << "[-:";
  } else {
    if (last_op_ != -1) {
      ofs << "[" << last_op_ << ":";
    } else {
      if (how_created_ == CREATE_MUTATE)
	ofs << "[m:";
      else if (how_created_ == CREATE_RECOMBINE)
	ofs << "[r:";
      else 
	ofs << "[?:";
    }
  }

  ofs << num_mutate_ << ":" << num_recombine_ << "]";
}


void IndHistory :: print_created(ostream& ofs)
{
  if (how_created_ == CREATE_RANDOM) {
    ofs << "[*(" << age_ << "):";
  } else if (how_created_ == CREATE_INIT) {
    ofs << "[-(" << age_ << "):";
  } else {
    if (how_created_ == CREATE_MUTATE)
      ofs << "[m(" << age_ << "):";
    else if (how_created_ == CREATE_RECOMBINE)
      ofs << "[r(" << age_ << "):";
    else 
      ofs << "[?(" << age_ << "):";
  }

  ofs << num_mutate_ << ":" << num_mutate_ << "]";
}


void IndHistory :: print_ops(ostream& ofs)
{
  if (op_count_.size() == 0)
    return;

  ofs << "[";
  for (vector<int>::size_type i=0; i<op_count_.size()-1; i++) {
    ofs << op_count_[i] << " ";
  }

  ofs << op_count_[op_count_.size()-1] << "]";
}


void IndHistory :: print_fitnesses(ostream& ofs)
{
  if (fitness_.size() == 1) {
    ofs << fitness_[0] << " ";
  } else {
    for (vector<double>::const_iterator iter = fitness_.begin();
	 iter != fitness_.end(); ++iter) {
      ofs << (*iter) << " ";
    }
  }
}

void IndHistory :: print(ostream& ofs)
{
  print_fitnesses(ofs);
  print_created(ofs);
}

void IndHistory :: print_f1(ostream& ofs)
{
  ofs << fitness_[0] << " ";
  print_created(ofs);
}


istream& IndHistory :: read(istream& istr)
{
  string line;
  vector<string> words;

  // Line 1: Read misc data
  getline(istr, line);
  split(line, ' ', words);
  if (words.size() != 6) {
    cout << "ind_history :: read() - error reading line 1:" << line << endl;
    while (1) ;
    throw runtime_error("ind_history :: read() - error reading line 1:" + line);
  }

  age_ = from_string<int>(words[0]);
  num_mutate_ = from_string<int>(words[1]);
  num_recombine_ = from_string<int>(words[2]);
  how_created_ = from_string<int>(words[3]);
  last_op_ = from_string<int>(words[4]);

  int fit_size = from_string<int>(words[5]);
  /*
  int fit_prev_size = from_string<int>(words[6]);
  int op_count_size = from_string<int>(words[7]);
  int complex_size = from_string<int>(words[8]);
  */

  fitness_.clear();
  fitness_.resize(0);
  getline(istr, line);
  split(line, ' ', words);
  for (int i=0; i<fit_size; i++) {
    double val = from_string<double>(words[i]);
    fitness_.push_back(val);
  }

  /*  
  fitness_prev_.clear();
  fitness_prev_.resize(0);
  for (int i=0; i<fit_prev_size; i++) {
    double val;
    istr >> val;
    fitness_prev_.push_back(val);
  }

  op_count_.clear();
  op_count_.resize(0);
  for (int i=0; i<op_count_size; i++) {
    int val;
    istr >> val;
    op_count_.push_back(val);
  }

  complexity_.clear();
  complexity_.resize(0);
  for (int i=0; i<complex_size; i++) {
    int val;
    istr >> val;
    complexity_.push_back(val);
  }
  */

  return istr;
}


ostream& IndHistory :: write(ostream& ostr) const
{
  ostr << age_ << " ";
  ostr << num_mutate_ << " ";
  ostr << num_recombine_ << " ";
  ostr << how_created_ << " ";
  ostr << last_op_ << " ";
  ostr << fitness_.size() << " ";
  /*
  ostr << fitness_prev_.size() << " ";
  ostr << op_count_.size() << " ";
  ostr << complexity_.size();
  */
  ostr << endl;

  {
    vector<double>::const_iterator iter = fitness_.begin();
    ostr << (*iter);
    for (++iter; iter != fitness_.end(); ++iter) {
      ostr << " " << (*iter);
    }
    ostr << endl;
  }

  /*
  {
    vector<double>::const_iterator iter = fitness_prev_.begin();
    ostr << (*iter);
    for (++iter; iter != fitness_prev_.end(); ++iter) {
      ostr << " " << (*iter);
    }
    ostr << endl;
  }

  if (!op_count_.empty()) {
    vector<int>::const_iterator iter = op_count_.begin();
    ostr << (*iter);
    for (++iter; iter != op_count_.end(); ++iter) {
      ostr << " " << (*iter);
    }
    ostr << endl;
  }

  if (!complexity_.empty()) {
    vector<double>::const_iterator iter = complexity_.begin();
    ostr << (*iter);
    for (++iter; iter != complexity_.end(); ++iter) {
      ostr << " " << (*iter);
    }
    ostr << endl;
  }
  */

  return ostr;
}

