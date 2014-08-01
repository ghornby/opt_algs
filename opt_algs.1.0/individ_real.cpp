/***

This file is part of the OptAlgs C++ library originally developed
by Gregory Hornby.

This code is released under the BSD-3 license:
  http://opensource.org/licenses/BSD-3-Clause
See the file "license.txt" in the root directory for full details.

***/

/*******************************************************

  File: individ_real.cpp

********************************************************/

#include "individ_real.h"

#include "cmn_logger.h"
#include "cmn_random.h"
#include "cmn_string_proc.h"


#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace std;
using namespace CmnClass;


#define OP_MUTATE_ALL 0
#define OP_MUTATE_SOME 1
#define OP_RECOMB_ALL 2
#define OP_RECOMB_SOME 3
#define OP_RECOMB_SWAP 4



#define INDIVID_REAL_VERSION 1

/***********************************************************/

const string Individ_Real :: class_name_("Individ_Real");


ostream& operator<<(ostream& os, const Individ_Real& a)
{
  return a.write(os);
}

istream& operator>>(istream& is, Individ_Real& a)
{
  return a.read(is);
}


/***********************************************************/

// Individ_real class functions.

Individ_Real :: Individ_Real(int num)
  : Individual(), num_genes_(num)
{
  gene_vec_.resize(num_genes_);
  max_val_.resize(num_genes_);
  min_val_.resize(num_genes_);
  init_max_.resize(num_genes_);
  init_min_.resize(num_genes_);
  set_minmax(0.0, 1.0);
  set_init_minmax(0.0, 1.0);
}


Individ_Real :: Individ_Real(const Individual *sample)
  : Individual()
{
  Individ_Real *s = (Individ_Real*)sample;

  num_genes_ = s->num_genes_;
  gene_vec_.resize(num_genes_);

  max_val_ = s->max_val_;
  min_val_ = s->min_val_;
  init_max_ = s->init_max_;
  init_min_ = s->init_min_;

  if (sample != 0) {
    num_genes_ = s->num_genes_;
  }
}


Individ_Real :: ~Individ_Real()
{
}


Individ_Real* Individ_Real :: new_instance() const
{
  //  CmnClass::Logger::Msg() << "individ_real :: new_instance() - " << num_genes_ << endl;
  Individ_Real* n = new Individ_Real(num_genes_);
  n->duplicate_settings(this);
  return n;
}


// Initialize variables.
void Individ_Real :: clear()
{
  Individual::clear();

  for (vector<double>::size_type i=0; i<gene_vec_.size(); i++) {
    gene_vec_[i] = init_min_[i];
  }
}


bool Individ_Real :: same(const Individ_Real *individ2) const
{
  int result = true;

  if (num_genes_ != individ2->num_genes_) {
    if (is_print_log_) {
      CmnClass::Logger::Msg() << "individual :: same() - number of genes different " << num_genes_
	   << " != " << individ2->num_genes_ << ".\n" << endl;
    }
    result = false;
  }

  if (!result) {
    while (1) ;
  }

  return result;
}


void Individ_Real :: set_minmax(vector<double>& min_val_s,
				vector<double>& max_val_s)
{
  min_val_ = min_val_s;
  max_val_ = max_val_s;
}


void Individ_Real :: set_minmax(double min, double max)
{
  for (vector<double>::iterator iter = min_val_.begin();
       iter != min_val_.end(); iter++) {
    *iter = min;
  }

  for (vector<double>::iterator iter = max_val_.begin();
       iter != max_val_.end(); iter++) {
    *iter = max;
  }
}


void Individ_Real :: set_init_minmax(double min, double max)
{
  for (vector<double>::iterator iter = init_min_.begin();
       iter != init_min_.end(); iter++) {
    *iter = min;
  }

  for (vector<double>::iterator iter = init_max_.begin();
       iter != init_max_.end(); iter++) {
    *iter = max;
  }
}


void Individ_Real :: set_init_minmax(vector<double>& mins,
				     vector<double>& maxs)
{
  init_min_ = mins;
  init_max_ = maxs;
}



void Individ_Real :: zero_genes()
{
  for (int i=0; i<num_genes_; i++) {
    gene_vec_[i] = 0;
  }
}  



int Individ_Real :: get_num_genes() const
{
  return num_genes_;
}


bool Individ_Real :: set_num_genes(int num)
{
  if (num <= 0) {
    return false;
  }

  num_genes_ = num;
  gene_vec_.resize(num);
  max_val_.resize(num_genes_);
  min_val_.resize(num_genes_);
  init_max_.resize(num_genes_);
  init_min_.resize(num_genes_);
  set_init_minmax(0.0, 1.0);
  set_minmax(0.0, 1.0);

  return true;
}

const vector<double>& Individ_Real :: get_genes() const
{
  return gene_vec_;
}

/**
 * Sets the genes for this individual. Genes are copied from source.
 *
 * @param genes Gene vector to be copied to this individual
 */
void Individ_Real::set_genes(std::vector<double>& genes) {
	gene_vec_ = genes;
}


// For this class, the "features" are exactly the genotype:
unsigned int Individ_Real :: num_features() const
{
  return gene_vec_.size();
}
vector<double> Individ_Real :: get_features() const
{
  return gene_vec_;
}
ostream& Individ_Real :: write_features(ostream& ostr) const
{
  CmnClass::write_vec(ostr, gene_vec_);
  return ostr;
}




void Individ_Real :: make_random()
{
  Individual::make_random();

  for (vector<double>::size_type i=0; i<gene_vec_.size(); i++) {
    gene_vec_[i] = random_double()*(init_max_[i]-init_min_[i]) + init_min_[i];
  }
}


void Individ_Real :: duplicate_settings(const Individ_Real *ind)
{
  Individual::duplicate_settings(ind);
  num_genes_ = ind->num_genes_;
  gene_vec_.resize(num_genes_);
  min_val_ = ind->min_val_;
  max_val_ = ind->max_val_;
  init_min_ = ind->init_min_;
  init_max_ = ind->init_max_;  
}
void Individ_Real :: duplicate_settings(const Individual *individ2)
{
  duplicate_settings((const Individ_Real*)individ2);
}



void Individ_Real :: duplicate(const Individ_Real *ind)
{
  Individual::duplicate(ind);
  gene_vec_ = ind->gene_vec_;
  min_val_ = ind->min_val_;
  max_val_ = ind->max_val_;
  init_min_ = ind->init_min_;
  init_max_ = ind->init_max_;
}
void Individ_Real :: duplicate(const Individual *ind)
{
  duplicate((const Individ_Real*)ind);
}



bool Individ_Real :: mutate(double scale)
{
  // This is only called by a failed recombination.

  //   double basic_scale = 0.05; // Original
  //  double basic_scale = 0.05 * scale; // 4.1 & 4.2
  //  double basic_scale = 0.04 * scale; // 4.3
  //  double basic_scale = 0.06 * scale; // 4.4
  //  double basic_scale = 0.05 * scale; // 4.5
  double basic_scale = 0.005 * scale * scale; // 4.5

  //  cout << "Individ_Real :: mutate(d)" << endl;


  Individual::mutate();

  int mutated = false;

  // Decide how many genes to mutate.

  /*
  //  const int k = 3; // Expers 1.1 - 4.1
  const int k = gene_vec_.size()*0.34; // Exper 4.2

  int max = (k <= (int)gene_vec_.size() ? k : (int)gene_vec_.size()); // Less if there are less genes.

  int num = 1; // Randomly select from [1-max]
  if (max >= 2) {
    //    num = random_int(max) + 1; // expers 1.1 - 4.1
    num = random_int(max) + 2; // Exper 4.2 ; note: perhaps allow a min of 1.
  }
  */

  int num = 0.2 * scale * random_norm() * gene_vec_.size(); // Exper 4.6
  num = abs(num) + 1;

  //  cout << "Num mutate: " << num << "   (scale:" << scale << ")" << endl;

  // Perform the mutations.
  for (int i=0; i<num; i++) {
    int index = 0;
    if (num_genes_ > 1) {
      index = random_int(num_genes_);
    }      
    double dist = max_val_[index] - min_val_[index];
    gene_vec_[index] += random_norm()*basic_scale*dist;
    if (gene_vec_[index] < min_val_[index]) {
      gene_vec_[index] = min_val_[index];
    } else if (gene_vec_[index] > max_val_[index]) {
      gene_vec_[index] = max_val_[index];
    }
    mutated = true;
  }


#ifdef ORIGINAL
  // experiments: 1.1 - 3.5
  // Decide how many genes to mutate.
  //  const int k = gene_vec_.size()*0.25; // At most k.  3?
  const int k = 3; // At most k.  3?
  int max = (k <= (int)gene_vec_.size() ? k : (int)gene_vec_.size()); // Less if there are less genes.
  int num = 1; // Randomly select from [1-max]
  if (max >= 2) {
    num = random_int(max) + 1;
  }

  // Perform the mutations.
  for (int i=0; i<num; i++) {
    int index = 0;
    if (num_genes_ > 1) {
      index = random_int(num_genes_);
    }      
    double dist = max_val_[index] - min_val_[index];
    gene_vec_[index] += random_norm()*basic_scale*dist;
    if (gene_vec_[index] < min_val_[index]) {
      gene_vec_[index] = min_val_[index];
    } else if (gene_vec_[index] > max_val_[index]) {
      gene_vec_[index] = max_val_[index];
    }
    mutated = true;
  }
#endif

  return mutated;
}

bool Individ_Real :: mutate()
{
  return mutate(1.0); // 0.1
  //  return mutate(0.1 * random_norm());
}



void Individ_Real :: get_unique(vector<Individual*>& inds2,
				vector<Individ_Real*>& indreals)
{
  indreals.resize(0);

  for (vector<Individual*>::size_type i=0; i<inds2.size(); i++) {
    Individ_Real* ind2 = (Individ_Real*)inds2[i];

    vector<double>::size_type j=0;
    for (; j<gene_vec_.size(); j++) {
      if (gene_vec_[j] != ind2->gene_vec_[j]) {
	break;
      }
    }
    if (j < gene_vec_.size()) {
      indreals.push_back(ind2);
    }
  }
}


bool Individ_Real :: mutate(double layer_scale, vector<Individual*>& inds2)
{
  // This is the main mutation function.

  if (inds2.size() == 0) {
    return mutate(layer_scale);
  }
  //  return mutate(layer_scale);

  /*
  cout << "Individ_Real :: mutate(d, vector<Ind*>&) :: layer_scale:" << layer_scale
       << ", #i2:" << inds2.size() << endl;
  */


  static unsigned int num_calls = 0;
  static unsigned int num_fails1 = 0;
  static unsigned int num_fails2 = 0;
  static unsigned int num_fails3 = 0;

  num_calls++;


  vector<Individ_Real*> p2s;
  get_unique(inds2, p2s);
  if (p2s.size() == 0) {
    //    CmnClass::Logger::Msg() << "Mutation :: error, all in tournament are the same.\n";
    num_fails1++;
    /*
    cout << "individ_real :: mutate(double, vec<Ind*>) - failed to mutate #1 ONE:"
	 << num_calls << "  ::  " << num_fails1 << " : " << num_fails2 << " : " << num_fails3 << endl;
    */
    return mutate(layer_scale);
  }


  // First figure out the range of values for each gene based on the
  // gene values in the tournament.
  bool is_all_same = true;
  vector<double> min = gene_vec_;
  vector<double> max = gene_vec_;
  for (vector<Individ_Real*>::size_type i=0; i<p2s.size(); i++) {
    for (vector<double>::size_type j=0; j<gene_vec_.size(); j++) {
      if (p2s[i]->gene_vec_[j] < min[j]) {
	min[j] = p2s[i]->gene_vec_[j];
	is_all_same = false;
      } else if (p2s[i]->gene_vec_[j] > max[j]) {
	max[j] = p2s[i]->gene_vec_[j];
	is_all_same = false;
      }
    }
  }

  if (is_all_same) {
    num_fails2++;
    /*
    cout << "individ_real :: mutate(double, vec<Ind*>) - failed to mutate #2 TWO:"
	 << num_calls << "  ::  " << num_fails1 << " : " << num_fails2 << " : " << num_fails3 << endl;
    */
    return mutate(layer_scale); // Exper 5.8
    //    return mutate(layer_scale/2.0);
  }


  unsigned int greater_range_ctr = 0;
  vector<unsigned int> mutatable_indexes;
  for (unsigned int i=0; i<gene_vec_.size(); i++) {
    if (min[i] < max[i]) {
      mutatable_indexes.push_back(i);
      double d = max[i] - min[i];
      if (d > 0.25 * (max_val_[i] - min_val_[i])) {
	greater_range_ctr++;
      }
    }
  }

  /*
  if (greater_range_ctr > 0) {
    cout << "Greater range ctr = " << greater_range_ctr << endl;
  }
  */

  // Randomize order, and the first index will ALWAYS be mutated:
  CmnClass::shuffle(mutatable_indexes);

  /*
  if (mutatable_indexes.size() < 10) {
    cout << "Mut indexes: ";
    for (vector<unsigned int>::size_type i=0; i<mutatable_indexes.size(); i++) {
      cout << "  " << mutatable_indexes[i];
    }
    cout << endl;
  }
  */

  bool is_mutated = false;
  unsigned int num_diff = 0;

  // Randomly determine how many genes/params to mutate:
  //  int num = 0.2 * layer_scale * random_norm() * gene_vec_.size(); // Exper 5.8
  /*
  int num = 0.2 * layer_scale * random_norm() * gene_vec_.size(); // Exper 5.8
  num = abs(num) + 1;
  if (num >= mutatable_indexes.size()) {
    num = mutatable_indexes.size();
  }
  */

    //    const double scale2 = 0.4; // Exper 5.10 - 5.11
    //    const double scale2 = 0.2; // Exper 5.12
    //    const double scale2 = 0.1; // Exper 5.13
    //    const double scale2 = 0.5; // Exper 5.14
    //    const double scale2 = 0.75; // Exper 5.15
    //    const double scale2 = 0.5; // Exper 5.16
    //    const double scale2 = 0.25; // Exper 6.13
    //    const double scale2 = 0.25; // Expers 5.3 and 6.14
    //    const double scale2 = 0.4; // Expers 5.3 and 6.14
    //    const double scale2 = 0.8; // Expers 5.3 and 6.14
  //    const double scale2 = 0.4; // Expers 5.3 and 6.14 and 7.1


  const double mutate_threshold = 0.7;
  const double scale2 = 1.0;
  const double max_range_perc = 1.0;

  for (vector<unsigned int>::size_type i=0; i<mutatable_indexes.size(); i++) {
  //  for (vector<unsigned int>::size_type i=0; i<num; i++) {
    if ((i != 0) && (random_double() < mutate_threshold)) { // 5.17
      //    if (random_double() < 0.5) { // 5.21
      //    if (random_double() < 0.4) { // 5.18
      //    if (random_double() < 0.75 * layer_scale) { // 5.19 
      continue;
    }

    // !!!!! Instead of adding to gene_vec_[i]; should this be from the
    // !!!!!   middle? or the average?  or should the alg keep track
    // !!!!    of range below and range above and scale differently in
    // !!!!!   each direction?


    unsigned int index = mutatable_indexes[i];
    double rand_dist = random_norm();


    // !!!!! Switch back to: (?)
    /*
    double range = max[index] - min[index];
    gene_vec_[index] += scale2 * layer_scale * rand_dist * range; // Exper: 5.3; 6.16
    */

    const double range_max = max_range_perc * (max_val_[index] - min_val_[index]);

    if (gene_vec_[index] == min_val_[index]) {
      double range = max[index] - gene_vec_[index];
      if (range > range_max) {
	range = range_max;
      }
      gene_vec_[index] += abs(scale2 * layer_scale * rand_dist * range);
      
    } else if (gene_vec_[index] == max_val_[index]) {
      double range = gene_vec_[index] - min[index];
      if (range > range_max) {
	range = range_max;
      }
      gene_vec_[index] -= abs(scale2 * layer_scale * rand_dist * range);

    } else {
      //      gene_vec_[index] += scale2 * scale * range * dist; // 5.10
      //	gene_vec_[index] += 0.5 * scale2 * scale * range * dist; // 5.11
      double range_2_max = max[index] - gene_vec_[index];
      double range_2_min = gene_vec_[index] - min[index];

      if (CmnClass::random_double() < 0.5) {
	double range = range_2_max; // Exper: 5.21
	if (range_2_max < range_2_min) {
	  range = (range_2_max + range_2_min) / 2.0;
	}
	if (range > range_max) {
	  range = range_max;
	}
	gene_vec_[index] += abs(scale2 * layer_scale * rand_dist * range); // 5.12
	
      } else {
	double range = range_2_min;  // Exper: 5.21
	if (range_2_min < range_2_max) {
	  range = (range_2_max + range_2_min) / 2.0;
	}
	if (range > range_max) {
	  range = range_max;
	}
	gene_vec_[index] -= abs(scale2 * layer_scale * rand_dist * range); // 5.12
      }
    }


    if (gene_vec_[index] < min_val_[index]) {
      gene_vec_[index] = min_val_[index];
    } else if (gene_vec_[index] > max_val_[index]) {
      gene_vec_[index] = max_val_[index];
    }
    is_mutated = true;
  }

  if (!is_mutated) {
    num_fails3++;
    /*
    cout << "individ_real :: mutate(double, vec<Ind*>) - failed to mutate #3 THREE:"
	 << num_calls << "  ::  " << num_fails1 << " : " << num_fails2 << " : " << num_fails3
	 << " (#diff:" << num_diff << ")" << endl;
    */
  }

  if (!is_mutated) {
    //    mutate(layer_scale/2.0);
    mutate(layer_scale);
  }


  Individual::mutate(layer_scale);

  return true;
}

bool Individ_Real :: mutate(vector<Individual*>& inds2)
{
  return mutate(1.0, inds2);
}





// Clips correctly.
void Individ_Real :: recombine_linear_clip(double dist, vector<double>& parent1,
					   vector<double>& parent2)
{
  if (dist < 0.0) {
    // Extrapolating past P1, so need to make sure we don't go out of bounds.
    for (vector<double>::size_type i=0; i<parent1.size(); i++) {
      double val = gene_vec_[i] + (parent2[i] - parent1[i]) * dist;
      if (val > max_val_[i]) {
	dist = (max_val_[i] - gene_vec_[i]) / (parent2[i] - parent1[i]);
      } else if (val < min_val_[i]) {
	dist = (min_val_[i] - gene_vec_[i]) / (parent2[i] - parent1[i]);
      }
    }

  } else if (dist > 1.0) {
    // Extrapolating past P2, so need to make sure we don't go out of bounds.
    for (vector<double>::size_type i=0; i<parent1.size(); i++) {
      double val = gene_vec_[i] + (parent2[i] - parent1[i]) * dist;
      if (val > max_val_[i]) {
	dist = (max_val_[i] - gene_vec_[i]) / (parent2[i] - parent1[i]);
      } else if (val < min_val_[i]) {
	dist = (min_val_[i] - gene_vec_[i]) / (parent2[i] - parent1[i]);
      }
    }
  }

  /*
  for (vector<double>::size_type i=0; i<parent1.size(); i++) {
    //    parent1[i] += (parent2[i] - parent1[i]) * dist;
    gene_vec_[i] += (parent2[i] - parent1[i]) * dist;
    if (gene_vec_[i] < min_val_[i] - 0.000000001) {
      CmnClass::Logger::Msg() << "gene_vec_[" << i << "] = " << gene_vec_[i]
	   << " less than min: " << min_val_[i] << endl;
      while (1) ;
    } else if (gene_vec_[i] > max_val_[i] + 0.0000000001) {
      CmnClass::Logger::Msg() << "gene_vec_[" << i << "] = " << gene_vec_[i]
	   << " greater than max: " << max_val_[i] << endl;
      while (1) ;
    }
  }
  */

  for (vector<double>::size_type i=0; i<gene_vec_.size(); i++) {
    gene_vec_[i] += (parent2[i] - parent1[i]) * dist;

    // Needed to keep gene values inside min-max limits since the above
    // equations can suffer from precision errors.
    if (gene_vec_[i] < min_val_[i]) {
      gene_vec_[i] = min_val_[i];
    } else if (gene_vec_[i] > max_val_[i]) {
      gene_vec_[i] = max_val_[i];
    }
  }
}


void Individ_Real :: recombine_linear(double dist, vector<double>& parent1,
				      vector<double>& parent2)
{
  // Values 0=>1 go from this/parent1 => parent2.
  // Values 1=>inf go beyond parent2 away from parent1
  // values 0=> -inf go beyond parent1 away from parent2
  for (vector<double>::size_type i=0; i<gene_vec_.size(); i++) {
    gene_vec_[i] += (parent2[i] - parent1[i]) * dist;
    if (gene_vec_[i] < min_val_[i]) {
      gene_vec_[i] = min_val_[i];
    } else if (gene_vec_[i] > max_val_[i]) {
      gene_vec_[i] = max_val_[i];
    }
  }
}


bool Individ_Real :: recombine_rand2(double scale, Individ_Real *parent2)
{
  if (num_genes_ < 1) {
    return false;
  }

  Individual::recombine(parent2);

  op_incr_used(OP_RECOMB_ALL);

  vector<double> tmp(gene_vec_);  
  double dist = scale * (2.0 * random_double() - 1.0);
  recombine_linear(dist, gene_vec_, parent2->gene_vec_);
  for (int i=0; i<num_genes_; i++) {
    if (random_double() < 0.25) {
      gene_vec_[i] = tmp[i];
    }
  }

  return true;
}
bool Individ_Real :: recombine_rand2(Individ_Real *parent2)
{
  return recombine_rand2(1.0, parent2);
}

bool Individ_Real :: recombine_rand2(double scale, Individual *parent2)
{
  //  return recombine((Individ_Real*)parent2);
  return recombine_rand2(scale, (Individ_Real*)parent2);
}
bool Individ_Real :: recombine_rand2(Individual *parent2)
{
  return recombine_rand2(1.0, (Individ_Real*)parent2);
}





bool Individ_Real :: recombine(double scale, Individ_Real *parent2)
{
  if (num_genes_ < 1) {
    return false;
  }

  //  Individual::recombine(parent2);
  op_incr_used(OP_RECOMB_ALL);    

  /*
  Individ_Real* ind2 = parent2;
  double msize = 0.75 * random_norm();
  for (int i=0; i<num_genes_; i++) {
    double dist = gene_vec_[i] - ind2->gene_vec_[i];
    gene_vec_[i] += msize*dist;
    if (gene_vec_[i] < min_val_[i]) {
      gene_vec_[i] = min_val_[i];
    } else if (gene_vec_[i] > max_val_[i]) {
      gene_vec_[i] = max_val_[i];
    }
  }
  */

  vector<double> old_val(gene_vec_);  

  Individ_Real* ind2 = parent2;

  //  double msize = 0.75 * random_norm(); // Ratio to perturb gene values.
  double msize = 0.5 * random_norm(); // Ratio to perturb gene values.
  msize *= scale;
  for (int i=0; i<num_genes_; i++) {
    double dist = gene_vec_[i] - ind2->gene_vec_[i]; // Distance btw P1 and P2.
    gene_vec_[i] += msize*dist; // New gene value.

    // Clip if it exceeds min/max:
    if (gene_vec_[i] < min_val_[i]) {
      gene_vec_[i] = min_val_[i];
    } else if (gene_vec_[i] > max_val_[i]) {
      gene_vec_[i] = max_val_[i];
    }

    /*
    // Restore to old value for some percentage of them:
    if (random_double() < 0.5) {
      gene_vec_[i] = old_val[i];
    }
    */
  }

  return true;
}
bool Individ_Real :: recombine(Individ_Real *parent2)
{
  return recombine(1.0, parent2);
}


bool Individ_Real :: recombine(double scale, Individual *parent2)
{
  return recombine(scale, (Individ_Real*)parent2);
}
bool Individ_Real :: recombine(Individual *parent2)
{
  return recombine(1.0, (Individ_Real*)parent2);
}



// This is the main recombination function.
bool Individ_Real :: recombine(double layer_scale, vector<Individual*>& inds2)
{
  if (num_genes_ < 1) {
    return false;
  }

  if (inds2.size() == 0) {
    return false;
  }


  //  cout << "individ_real :: recombine_rand(layer_scale, vector<Ind>) - called" << endl;

  //  while (1) ;


  //  cout << "individ_real :: recombine_rand(vector<Ind>) - layer_scale: " << layer_scale << endl;

  //  Individual::recombine(inds2[0]);
  return recombine(inds2[0]);
  op_incr_used(OP_RECOMB_ALL);    

  vector<Individ_Real*> p2s;
  get_unique(inds2, p2s);
  if (p2s.size() == 0) {
    //    CmnClass::Logger::Msg() << "Recombine :: error, all in tournament are the same.\n";
    return mutate();
  }

  int recomb_type = 1;
  /*
  if (random_double() < 0.5) {
    recomb_type = 6;
  }
  */


  if (recomb_type == 1) {
    // Recomb 1
    Individ_Real* ind2 = (Individ_Real*)p2s[0];
    if (p2s.size() > 1) {
      int index = random_int(p2s.size());
      ind2 = (Individ_Real*)p2s[index];
    }



    //    double dist = random_norm();
    //    double dist = 0.25 * random_norm();
    // 2.3?:   double dist = random_double() * 2.0 - 1.0;
    // 2.4:   double dist = random_double() * 2.25 - 0.75;
    /* 2.5:
    double dist = random_double() * 1.8 - 0.9;
    if (dist > 0) {
      dist *= 1.25/0.9;
    }
    */
    /* 2.7:
    double dist = random_double();
    if (random_double() < 0.5) {
      dist = -abs(random_norm());
    }
    */
    /*
    do {
      dist = 0.75 * random_norm();
    } while (dist > 1.0);
    */

    // 2.8:
    //    double dist = (random_double() * 2.0 - 1.0 + random_norm()) / 2.0;

    // 2.9:
    //    double dist = (random_double() * 2.0 - 1.0 + 1.5 * random_norm()) / 2.0;

    // 2.10:
    //    double dist = (random_double() * 2.0 - 1.0 + 0.8 * random_norm()) / 2.0;

    // 2.11:
    //    double dist = (random_double() * 2.0 - 1.0 + 1.2 * random_norm()) / 2.0;

    // 2.12:
    //    double dist = (random_double() * 2.0 - 1.0 + 0.9 * random_norm()) / 2.0;

    // 2.13:(age-gap:3) & 2.14:(age-gap:2); and 2.15-2.17
    //    double dist = (random_double() * 2.0 - 1.0 + 0.7 * random_norm()) / 2.0;
    //    dist *= layer_scale;


    //    const double weight_uniform = 0.5; // Exper 3.1 - 5.21
    //    const double weight_uniform = 0.25; // Exper 6.1
    //    const double weight_uniform = 0.75; // Exper 6.2
    //    const double weight_uniform = 0.4; // Exper 6.3-6.6
    //    const double scale_norm = 0.75; // Exper 6.3
    //    const double scale_norm = 1.25; // Exper 6.4
    //    const double scale_norm = 1.4; // Exper 6.5
    //    const double scale_norm = 0.8; // Exper 6.6

    //    const double weight_uniform = 0.6; // Exper 6.7
    //    const double scale_norm = 1.2; // Exper 6.7

    //    const double weight_uniform = 0.4; // Exper 6.11
    //    const double scale_norm = 0.5; // Exper 6.11

    //    const double weight_uniform = 0.0; // Exper 6.12
    //    const double scale_norm = 0.6; // Exper 6.12

    const double weight_uniform = 0.4; // Exper 7.X
    const double scale_norm = 0.9; // Exper 7.X


    double dist = weight_uniform * (random_double() * 2.0 - 1.0);
    dist += layer_scale * scale_norm * (1.0 - weight_uniform) * random_norm();
    //    dist *= layer_scale; // Exper 3.x? - 6.6

    //    cout << "individ_real :: recombine(layer_scale, vec<Ind*>) - dist:" << dist << endl;


    /* // 3.1
    // uniform part of distance:
    double d_uniform = random_double();  // Default is interpolating
    if (random_double() < 0.5) {
      m_uniform *= (-1.25); // Extrapolating K% of dist between P1 and P2
    }

    double d_norm = random_norm();

    double dist = 0.5*d_uniform + 0.5*d_norm;
    */

      
    //    double msize = random_norm();
    recombine_linear(dist, gene_vec_, ind2->gene_vec_);
    //    recombine_linear_clip(dist, gene_vec_, ind2->gene_vec_);  // <= doesn't seem to work as well as recombine_linear()
    /*
    for (int i=0; i<num_genes_; i++) {
      double dist = ind1->gene_vec_[i] - ind2->gene_vec_[i];
      gene_vec_[i] += msize*dist;
      if (gene_vec_[i] < min_val_[i]) {
	gene_vec_[i] = min_val_[i];
      } else if (gene_vec_[i] > max_val_[i]) {
	gene_vec_[i] = max_val_[i];
      }
    }
    */


  } else if (recomb_type == 2) {
    // Recomb 2
    vector<double> mut_range(gene_vec_.size());
    {
      Individ_Real* ind2 = (Individ_Real*)p2s[0];
      if (p2s.size() > 1) {
	ind2 = (Individ_Real*)p2s[random_int(p2s.size())];
      }

      for (vector<double>::size_type i=0; i<gene_vec_.size(); i++) {
	mut_range[i] = gene_vec_[i] - ind2->gene_vec_[i];
      }
    }


    Individ_Real* ind2 = (Individ_Real*)p2s[0];
    if (p2s.size() > 1) {
      ind2 = (Individ_Real*)p2s[random_int(p2s.size())];
    }

    //    double msize = random_double() * 2.0 - 1.0;
    //    double msize = random_double() * 1.8 - 0.9;
    double msize = layer_scale * 0.75 * random_norm();
    recombine_linear(msize, gene_vec_, ind2->gene_vec_);

    for (int i=0; i<num_genes_; i++) {
      gene_vec_[i] += layer_scale * 0.01 * random_norm() * random_norm() * mut_range[i];
      if (gene_vec_[i] < min_val_[i]) {
	gene_vec_[i] = min_val_[i];
      } else if (gene_vec_[i] > max_val_[i]) {
	gene_vec_[i] = max_val_[i];
      }
    }


  } else if (recomb_type == 3) {
    // Recomb 3
    Individ_Real* ind2 = (Individ_Real*)p2s[0];
    if (p2s.size() > 1) {
      ind2 = (Individ_Real*)p2s[random_int(p2s.size())];
    }

    //    double msize = random_double() * 2.0 - 1.0;
    //    double msize = random_double() * 1.8 - 0.9;
    double msize = layer_scale * 1.0 * random_norm();
    recombine_linear(msize, gene_vec_, ind2->gene_vec_);
    /*
    for (int i=0; i<num_genes_; i++) {
      double dist = gene_vec_[i] - ind2->gene_vec_[i];
      gene_vec_[i] += msize*dist;
      if (gene_vec_[i] < min_val_[i]) {
	gene_vec_[i] = min_val_[i];
      } else if (gene_vec_[i] > max_val_[i]) {
	gene_vec_[i] = max_val_[i];
      }
    }
    */

  } else if (recomb_type == 4) {
    // Recomb 4
    // First: randomly select two DIFFERENT parents:
    Individ_Real* ind1 = this;
    Individ_Real* ind2 = (Individ_Real*)p2s[0];
    if (p2s.size() > 1) {
      int index = random_int(p2s.size()+1);
      if (index != (int)p2s.size()) {
	ind1 = (Individ_Real*)p2s[index];
      } else {
	// ind1 = this. // Set above.
      }

      index = random_int(p2s.size());
      ind2 = (Individ_Real*)p2s[index];
      if (ind2 == ind1) {
	ind1 = this;
      }
    }

    //    double msize = 0.25 * random_norm();
    /*
    double msize = random_double() * 2.0 - 1.0;
    for (int i=0; i<num_genes_; i++) {
      double dist = ind1->gene_vec_[i] - ind2->gene_vec_[i];
      gene_vec_[i] += msize*dist;
      if (gene_vec_[i] < min_val_[i]) {
	gene_vec_[i] = min_val_[i];
      } else if (gene_vec_[i] > max_val_[i]) {
	gene_vec_[i] = max_val_[i];
      }
    }
    */
    double msize = random_double() * 2.0 - 1.0;
    msize *= layer_scale;
    recombine_linear(msize, ind1->gene_vec_, ind2->gene_vec_);


  } else if (recomb_type == 5) {
    // Recomb 5
    vector<double> mut_range(gene_vec_.size());
    {
      Individ_Real* ind2 = (Individ_Real*)p2s[0];
      if (p2s.size() > 1) {
	ind2 = (Individ_Real*)p2s[random_int(p2s.size())];
      }

      for (vector<double>::size_type i=0; i<gene_vec_.size(); i++) {
	mut_range[i] = gene_vec_[i] - ind2->gene_vec_[i];
      }
    }


    vector<double> gradient(gene_vec_.size());
    for (vector<double>::size_type j=0; j<gradient.size(); j++) {
      gradient[j] = gene_vec_[j];
    }
    int ctr = p2s.size();
    double num_pos = 1.0;
    for (vector<Individual*>::size_type i=0; i<p2s.size(); i++) {
      ctr -= 2;
      if (ctr == 0) {
	continue;
      } else if (ctr > 1) {
	num_pos += 1.0;
      }

      double mult = 1.0;
      if (ctr < 0.0) {
	mult = -1.0;
      }
      Individ_Real* ind2 = (Individ_Real*)p2s[i];
      for (vector<double>::size_type j=0; j<gradient.size(); j++) {
	gradient[j] += ind2->gene_vec_[j] * mult;
      }
    }


    //    double div = p2s.size() * (p2s.size() + 1.0) / 2.0;
    //    double div = num_pos;
    /*
    for (vector<double>::size_type j=0; j<gradient.size(); j++) {
      gradient[j] /= div;
    }
    */

    Individ_Real* ind1 = this;
    /*
    int index = random_int(inds2.size()+1);
    if (index < inds2.size()) {
      ind1 = (Individ_Real*)inds2[index];
    }
    */

    double dist = layer_scale * 0.75 * random_norm();
    //    double dist = 2.0*random_double() - 1.0;
    for (int i=0; i<num_genes_; i++) {
      gene_vec_[i] = ind1->gene_vec_[i] + gradient[i]*dist;
      if (gene_vec_[i] < min_val_[i]) {
	gene_vec_[i] = min_val_[i];
      } else if (gene_vec_[i] > max_val_[i]) {
	gene_vec_[i] = max_val_[i];
      }
    }

    for (int i=0; i<num_genes_; i++) {
      gene_vec_[i] += 0.01 * random_norm() * random_norm() * mut_range[i];
      if (gene_vec_[i] < min_val_[i]) {
	gene_vec_[i] = min_val_[i];
      } else if (gene_vec_[i] > max_val_[i]) {
	gene_vec_[i] = max_val_[i];
      }
    }

  } else if (recomb_type == 6) {
    // Recomb 6
    vector<double> mut_range(gene_vec_.size());
    {
      Individ_Real* ind2 = (Individ_Real*)p2s[0];
      if (p2s.size() > 1) {
	ind2 = (Individ_Real*)p2s[random_int(p2s.size())];
      }

      for (vector<double>::size_type i=0; i<gene_vec_.size(); i++) {
	mut_range[i] = gene_vec_[i] - ind2->gene_vec_[i];
      }
    }


    double mult = p2s.size();
    vector<double> gradient(gene_vec_.size());
    for (vector<double>::size_type j=0; j<gradient.size(); j++) {
      gradient[j] = gene_vec_[j] * mult;
    }
    for (vector<Individual*>::size_type i=0; i<p2s.size(); i++) {
      mult -= 2.0;
      if (mult == 0.0) {
	continue;
      }

      Individ_Real* ind2 = (Individ_Real*)p2s[i];
      for (vector<double>::size_type j=0; j<gradient.size(); j++) {
	gradient[j] += ind2->gene_vec_[j] * (double)mult;
      }
    }
    double div = p2s.size() * (p2s.size() + 1.0) / 2.0;
    for (vector<double>::size_type j=0; j<gradient.size(); j++) {
      gradient[j] /= div;
      //      gradient[j] = gradient[j] * 2.0 / div;
    }


    /*
    for (vector<double>::size_type i=0; i<gene_vec_.size();
	 i++) {
      CmnClass::Logger::Msg() << gradient[i] << ":" << gene_vec_[i] - p2s[0]->gene_vec_[i]
	   << ", ";
    }
    CmnClass::Logger::Msg() << endl;
    */

    /*
    Individ_Real* ind1 = this;
    int index = random_int(inds2.size()+1);
    if (index < inds2.size()) {
      ind1 = (Individ_Real*)inds2[index];
    }
    */

    double dist = layer_scale * 0.75 * random_norm();
    //    double dist = random_norm();
    //    double dist = random_double()*2.0 - 1.0;

    for (int i=0; i<num_genes_; i++) {
      gene_vec_[i] += dist * gradient[i];
      if (gene_vec_[i] < min_val_[i]) {
	gene_vec_[i] = min_val_[i];
      } else if (gene_vec_[i] > max_val_[i]) {
	gene_vec_[i] = max_val_[i];
      }
    }

    /*
    for (int i=0; i<num_genes_; i++) {
      gene_vec_[i] += 0.1 * random_norm() * random_norm() * mut_range[i];
      if (gene_vec_[i] < min_val_[i]) {
	gene_vec_[i] = min_val_[i];
      } else if (gene_vec_[i] > max_val_[i]) {
	gene_vec_[i] = max_val_[i];
      }
    }
    */

    // Decide how many genes to mutate.
    const int k = 5;
    int max = (k <= (int)gene_vec_.size() ? k : (int)gene_vec_.size());
    int num = 1;
    if (max > 2) {
      num = random_int(max-1) + 1;
    }    
    for (int i=0; i<num; i++) {
      int index = random_int(num_genes_);
      gene_vec_[index] += layer_scale * 0.02 * random_norm() * mut_range[index];
      if (gene_vec_[index] < min_val_[index]) {
	gene_vec_[index] = min_val_[index];
      } else if (gene_vec_[index] > max_val_[index]) {
	gene_vec_[index] = max_val_[index];
      }
    }
  }


  return true;
}
bool Individ_Real :: recombine(vector<Individual*>& inds2)
{
  return recombine(1.0, inds2);
}



int Individ_Real :: phenotype_distance(Individ_Real* ind2)
{
  int num = min(num_genes_, ind2->num_genes_);

  int dist = max(num_genes_, ind2->num_genes_) - num;
  for (int i=0; i<num; i++) {
    if (gene_vec_[i] != ind2->gene_vec_[i]) {
      num++;
    }
  }

  return dist;
}
int Individ_Real :: phenotype_distance(Individual* ind2)
{
  return phenotype_distance((Individ_Real*)ind2);
}


  //bool Individ_Real :: make_phenotype(void* arg)
bool Individ_Real :: make_phenotype(void*)
{
  return true;
}

bool Individ_Real :: make_phenotype()
{
  return true;
}


double Individ_Real :: distance_apart(const Individ_Real* ind2) const
{
  double sum = 0;

  for (vector<double>::size_type i=0; i<gene_vec_.size(); i++) {
    double diff = gene_vec_[i] - ind2->gene_vec_[i];
    sum += diff * diff;
  }

  return sqrt(sum);
}
double Individ_Real :: distance_apart(const Individual* ind2) const
{
  return distance_apart((const Individ_Real*)ind2);
}




istream& Individ_Real :: read(istream& istr)
{
  string line;
  vector<string> words;

  // 0: Read type and version number:
  getline(istr, line);
  split(line, ' ', words);
  if (words.size() != 2) {
    CmnClass::Logger::Msg() << "IndReal :: read() - error reading version:" << line << endl;

    string line2;
    getline(istr, line2);
    CmnClass::Logger::Msg() << "IndReal :: read() - Next line:" << line2 << endl;

    while (1) ;
    throw runtime_error("IndReal :: read() - error reading version:" + line);
  }

  if (words[0] != "IndRealV") {
    CmnClass::Logger::Msg() << "IndReal :: read() - error reading version[0]:" << line << endl;
    while (1) ;
    throw runtime_error("IndReal :: read() - error, expecting IndV, got:" + line);
  }

  int version = from_string<int>(words[1]);
  if (version != INDIVID_REAL_VERSION) {
    CmnClass::Logger::Msg() << "IndReal :: read() - error reading version[1]:" << line << endl;
    while (1) ;
    throw runtime_error("IndReal :: read() - error, bad version:" + words[1]);
  }  


  Individual::read(istr);

  read_vec(istr, gene_vec_, "num genes");
  num_genes_ = gene_vec_.size();

  read_vec(istr, max_val_, "max val");
  read_vec(istr, min_val_, "min val");
  read_vec(istr, init_max_, "init max");
  read_vec(istr, init_min_, "init min");

  return istr;
}


bool Individ_Real :: read(const char *fname)
{
  CmnClass::Logger::Msg() << "Individ_Real :: read() " << fname << endl;
  ifstream file(fname);
  try {
    read(file);
    //  } catch (AlpsError e) {
  } catch(...) {
    return false;
  }

  file.close();
  return true;
}


ostream& Individ_Real :: write(std::ostream& ostr) const
{
  ostr << "IndRealV " << INDIVID_REAL_VERSION << endl;

  Individual::write(ostr);

  int p = ostr.precision();
  ostr.precision(10);

  write_vec(ostr, gene_vec_);
  write_vec(ostr, max_val_);
  write_vec(ostr, min_val_);
  write_vec(ostr, init_max_);
  write_vec(ostr, init_min_);

  ostr.precision(p);
  return ostr;
}


bool Individ_Real :: write(const char *fname) const
{
  ofstream file(fname);
  write(file);
  return true;
}


/********************************************************/


