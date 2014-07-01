/***

This file is part of the OptAlgs C++ library originally developed
by Gregory Hornby.

This code is released under the BSD-3 license:
  http://opensource.org/licenses/BSD-3-Clause
See the file "license.txt" in the root directory for full details.

***/

/*******************************************************

  File: individual.cpp

********************************************************/


#include "individual.h"


#include "cmn_logger.h"
#include "cmn_random.h"
#include "cmn_string_proc.h"



#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <stdexcept>

using namespace CmnClass;
using namespace std;


const int INDIVIDUAL_VERSION = 1;

const int NUM_VARIATION_OPS = 10;


/***********************************************************/

ostream& operator<<(ostream& os, const Individual& a)
{
  return a.write(os);
}

istream& operator>>(istream& is, Individual& a)
{
  return a.read(is);
}

/***********************************************************/


int gNum_Fitness = 5;
int gNum_Operators = NUM_VARIATION_OPS;
double gRec_delta = 1.0;
double gMutate_size = 1.0;


int gNum_Mutations = 0;
int gNum_Mutations_Equal = 0;
int gNum_Mutations_Better = 0;
int gNum_Recombinations = 0;
int gNum_Recombinations_Equal = 0;
int gNum_Recombinations_Better = 0;


int gOp_Used[NUM_HISTORY_OPS];
int gOp_Used_Hist[NUM_HISTORY_OPS];

const int OP_STACK_SIZE = 20;
int gOp_Stack[OP_STACK_SIZE];
int gOp_Stack_Ptr = 0;


/***********************************************************/

const int Age_Take_Older = 1;
const int Age_Take_Average = 2;
const int Age_Take_Younger = 3;


/***********************************************************/


bool op_incr_used(int op_num)
{
  if (gOp_Stack_Ptr >= OP_STACK_SIZE)
    return false;

  gOp_Stack[gOp_Stack_Ptr++] = op_num;
  return true;
}

int op_get_last()
{
  if (gOp_Stack_Ptr > 0)
    return gOp_Stack[gOp_Stack_Ptr-1];
  return -1;
}


void op_reset(void)
{
  gOp_Stack_Ptr = 0;
}

void op_update(vector<int>& operation)
{
  int i;

  //  cout << "individual :: Op_update() - #" << Op_stack_ptr << endl;
  for (i=0; i<gOp_Stack_Ptr; i++) {
    //    cout << "individual :: op_update(): " << Op_stack[i] << endl;
    if (gOp_Stack[i] >= (int)operation.size()) {
      operation.resize(gOp_Stack[i]+1);
    }
    operation[gOp_Stack[i]]++;
  }
  /*
  if (better) {
    for (i=0; i<Op_stack_ptr; i++)
      Op_better[Op_stack[i]]++;
  }
  */
  gOp_Stack_Ptr = 0;
}


/***********************************************************/


inline double reduce(double val)
{
  return (floor((val) * 1.0E+8 + 0.5) / 1.0E+8);
}


double Rec_delta = 1.0;


double recombine_double(double val1, double val2)
{
  double diff;

  diff = Rec_delta * (val1 - val2);
  return reduce(val1 + diff * (CmnClass::random_double() * 2.0 - 1.0));
}


double recombine(double val1, double val2)
{
  if (val1 == 0.0) {
    if (random_int(2))
      return val1;
    else
      return val2;
  } else if (val2 == 0.0) {
    if (random_int(2))
      return val1;
    else
      return val2;
  }

  double diff;

  diff = Rec_delta * (val1 - val2);
  return reduce(val1 + diff * (random_double() * 2.0f - 1.0f));
}



void recombine_elliptic(vector<double>& parent1, vector<double>& parent2,
			vector<double>& child1)
{
  cerr << "alps_utils :: recombine_elliptic() - not implemented.\n";
  while (1) {
    ;
  }  

   /* child1 - the random vector in the circle */
   for (vector<double>::size_type i=0; i<child1.size(); i++) {
     child1[i] = random_norm();
   }


   /* normalize child1 to be a unit vector */
   double val = 0.0;
   for (vector<double>::size_type i=0; i<child1.size(); i++) {
     val += child1[i] * child1[i];
   }
   val = sqrt(val);

   for (vector<double>::size_type i=0; i<child1.size(); i++) {
     child1[i] /= val;
   }

   /* give it a random length */
   val = pow(random_double(), 1.0/(double)child1.size());
   for (vector<double>::size_type i=0; i<child1.size(); i++)
     child1[i] *= val;

   cout << "C1: ";
   for (vector<double>::size_type i=0; i<child1.size(); i++) {
     cout << child1[i] << " ";
   }  
   cout << "\n";

   
   /* Set up I, a and b, and d */
   
   /* d - the vector from p1 to p2 */
   //   Individual d = pop_string_temp(string_pop(parent1));
   //   Individual uu = pop_string_temp2(string_pop(parent1));

   vector<double> d(parent2);
   for (vector<double>::size_type i=0; i<parent1.size(); i++) {
     d[i] -= parent1[i];
   }

   val = 0.0;
   for (vector<double>::size_type i=0; i<parent1.size(); i++) {
     val += d[i] * d[i];
   }
   double I = sqrt(val);

   //   double a = I * 1.0 / sqrt((double)child1.size());
   double a = I * 1.0;
   double b = I * 0.3 / sqrt((double)child1.size());

   /* uu is projection of child1/u onto d */
   double sum = 0.0;
   for (vector<double>::size_type i=0; i<child1.size(); i++) {
     sum += child1[i] * d[i];
   }
	 
   val = sum / val;
   vector<double> uu(d);
   for (vector<double>::size_type i=0; i<parent1.size(); i++) {
     uu[i] *= val;
   }


   /* Now stretch child1 to proper shape and translate it */
   for (vector<double>::size_type i=0; i<child1.size(); i++) {
     child1[i] = b * (child1[i] - uu[i])
       + a * (uu[i] + parent1[i]);
   }
   
   
   cout << "P1: ";
   for (vector<double>::size_type i=0; i<parent1.size(); i++) {
     cout << parent1[i] << " ";
   }  
   cout << "\n";

   cout << "P2: ";
   for (vector<double>::size_type i=0; i<parent2.size(); i++) {
     cout << parent2[i] << " ";
   }  
   cout << "\n";



   cout << "d:  ";
   for (vector<double>::size_type i=0; i<d.size(); i++) {
     cout << d[i] << " ";
   }  
   cout << "\n";

   cout << "UU: ";
   for (vector<double>::size_type i=0; i<uu.size(); i++) {
     cout << setw(7) << uu[i] << " ";
   }  
   cout << "\n";



   cout << "C1: ";
   for (vector<double>::size_type i=0; i<child1.size(); i++) {
     cout << child1[i] << " ";
   }  
   cout << "\n";
}




const double INTEGER_RANGE = 100.0;
//const int INT_RANGE = 100;


int recombine_int(int val1, int val2)
{
  
  if (val1 == val2)
    return val1;

  if (val1 < val2) {
    int diff = val2 - val1;
    diff = random_int(2 * diff + 1);
//    printf("general :: recombine_int() - a. %d:%d -> diff=%d -> %d => ",
//    	   val1, val2, (val2-val1), diff);

    val2 -= diff;
    //    printf("%d.\n", val2);
    return val2;
  }

  int diff = val1 - val2;
  diff = random_int(2 * diff + 1);
  //  printf("general :: recombine_int() - b. %d:%d -> diff=%d -> %d => ",
  //	 val1, val2, (val1-val2), diff);
  
  val2 += diff;
  //  printf("%d,\n", val1);
  return val2;
}

double recombine_double_as_int(double val1, double val2)
{
  
  if (val1 == val2)
    return val1;

  int diff = int(val1 - val2);
  diff = random_int(2 * diff + 1) - diff;
  //  diff = rand()%(2 * diff + 1) - diff;

  val1 += double(diff);
  if (val1 > INTEGER_RANGE)
    val1 = INTEGER_RANGE;
  else if (val1 < 0.0)
    val1 = 0.0;
  return val1;
}


/***********************************************************/


bool Individual :: is_print_log_ = true;
const string Individual :: class_name_("Individual");


// Individual class functions.

Individual :: Individual()
  : is_deleted_(false), design_id_(0)
{
}


Individual :: Individual(const Individual *)
  : is_deleted_(false), design_id_(0)
{
  if (is_print_log_) {
    CmnClass::Logger::Msg() << "individual :: individual(individual* sample)" << endl;
  }
}


Individual :: ~Individual()
{
  if (is_deleted_) {
    CmnClass::Logger::Msg() << "individual :: delete() - warning, already deleted." << endl;
    while (1) ;
  }
  //  CmnClass::Logger::Msg() << " deleting individual " << design_id_ << endl;

  is_deleted_ = true;
}


bool Individual :: is_deleted() const
{
  return is_deleted_;
}

//string Individual :: get_class_name() const
//{
//  CmnClass::Logger::Msg() << "individual :: get_classname() - error, should be instantiated in subclass." << endl;
//  while (1) ;
//  return class_name_;
//}
//
//
//Individual* Individual :: new_instance() const
//{
//  return new Individual();
//}


// Clear/Initialize variables.
void Individual :: clear()
{
  history_self_.clear();
  history_parent1_.clear();
  history_parent2_.clear();

  evaluations_ = 0;
  design_id_ = 0;

  for (vector<double>::size_type i=0; i<features_.size(); i++) {
    features_[i] = 0;
  }
}


bool Individual :: valid()
{
  return true;
}


// fitness_hist_diff() - returns the different in fitness
// between this individual and its parent: (this fitness - parent fitness).
// parent is taken as the first parent.
double Individual :: fitness_hist_diff(int index)
{
  if ((index < 0) || (index >= gNum_Fitness)) {
    CmnClass::Logger::Msg() << "individual :: fitness_hist_diff() - error on fitness index: "
	 << index << "(0 to " << gNum_Fitness << ").\n";
    while (1) ;
  }

  return history_self_.get_fitness(index)
    - history_parent1_.get_fitness(index);
}

double Individual :: fitness_hist_diff()
{
  return (history_self_.get_fitness() - history_parent1_.get_fitness());
}



// Returns 1 if better fitness than parent(s)
// returns 0 if equal/inbetween parent(s)
// returns -1 if strictly less than parent(s).
int Individual :: better(bool is_maximizing)
{
  if (is_maximizing) {
    if ((history_self_.get_fitness() > history_parent1_.get_fitness()) &&
        (history_self_.get_fitness() > history_parent2_.get_fitness()))
      return 1;
    else if ((history_self_.get_fitness() == history_parent1_.get_fitness()) ||
             (history_self_.get_fitness() == history_parent2_.get_fitness()))
      return 0;
    else if ((history_self_.get_fitness() > history_parent1_.get_fitness()) ||
             (history_self_.get_fitness() > history_parent2_.get_fitness()))
      return 0;
    else
      return -1;
  }

    if ((history_self_.get_fitness() < history_parent1_.get_fitness()) &&
        (history_self_.get_fitness() < history_parent2_.get_fitness()))
      return 1;
    else if ((history_self_.get_fitness() == history_parent1_.get_fitness()) ||
             (history_self_.get_fitness() == history_parent2_.get_fitness()))
      return 0;
    else if ((history_self_.get_fitness() < history_parent1_.get_fitness()) ||
             (history_self_.get_fitness() < history_parent2_.get_fitness()))
      return 0;
    else
      return -1;
}


// Conditions for keeping this individual versus re-applying
// variation to the parent(s).
bool Individual :: keep(bool is_maximizing)
{
  if (is_maximizing) {
    if (history_self_.get_how_created() == CREATE_RANDOM) {
      if (history_self_.get_fitness() > 0.0) {
	/*
	  printf("individual :: keep() - init fitness %f > 0.0, keep.\n",
	  history_self_.get_fitness());
	*/
	return true;
      }
      /*
	printf("individual :: keep() - init fitness %f <= 0.0, don't keep.\n",
	history_self_.get_fitness());
      */
      return false;
    }
    
    if (history_self_.get_fitness() >= history_parent1_.get_fitness()*0.1) {
      /*
	printf("individual :: keep() - fitness %f > 0.1* %f keep.  (Created=%d).\n",
	history_self_.get_fitness(), history_parent1_.get_fitness(),
	history_self_.get_how_created());
      */
      return true;
    }

  } else {
    if (history_self_.get_how_created() == CREATE_RANDOM) {
      if (history_self_.get_fitness() < Worst_Fitness) {
	/*
	  printf("individual :: keep() - init fitness %f > 0.0, keep.\n",
	  history_self_.get_fitness());
	*/
	return true;
      }
      /*
	printf("individual :: keep() - init fitness %f <= 0.0, don't keep.\n",
	history_self_.get_fitness());
      */
      return false;
    }
    
    if (history_self_.get_fitness() <= history_parent1_.get_fitness()*100.0) {
      /*
	printf("individual :: keep() - fitness %f > 0.1* %f keep.  (Created=%d).\n",
	history_self_.get_fitness(), history_parent1_.get_fitness(),
	history_self_.get_how_created());
      */
      return true;
    }
  }


  /*
  printf("individual :: keep() - fitness %f < 0.1x %f, don't keep.\n",
	 history_self_.get_fitness(), history_parent1_.get_fitness());
  */

  return false;
}



// Returns 1 if this is better fitness than individ2
// returns 0 if this is equal to individ2
// returns -1 if this is strictly less than individ2.
int Individual :: compare_fitness(bool is_maximizing, Individual *individ2)
{
  //  double *fitness2 = individ2->history_self_.fitness;

  if (is_maximizing) {
    if (history_self_.get_fitness() > 
	individ2->history_self_.get_fitness()) {
      return 1; // This individual is better.
    }

  } else {
    if (history_self_.get_fitness() < 
	individ2->history_self_.get_fitness()) {
      return 1; // This individual is better.
    }
  }

  if (history_self_.get_fitness() == 
      individ2->history_self_.get_fitness()) {
    return 0; // Both are equal.
  }

  return -1; // The other individual is better.
}


//bool Individual :: same(Individual *individ2)
bool Individual :: same(const Individual *)
{
  // To be implemented by subclass.
  return true;
}


bool Individual :: created_by_variation()
{
  if ((history_self_.get_how_created() != CREATE_MUTATE) &&
      (history_self_.get_how_created() != CREATE_RECOMBINE)) {
    return false;
  }

  return (evaluations_ == 1);
}


double Individual :: get_fitness() const
{
  return history_self_.get_fitness();
}


void Individual :: get_fitness(vector<double>& fitness_vec) const
{
  fitness_vec = history_self_.get_fitness_vec();
}

const vector<double>& Individual :: get_fitness_vec() const
{
  return history_self_.get_fitness_vec();
}



void Individual :: set_fitness(bool maximize, const vector<double>& fitness_vec)
{
  int result;

  /*
  CmnClass::Logger::Msg() << "individual() :: set_fit() " << fitness_vec[0]
       << " evals=" << evaluations_ << endl;
  */

  history_self_.backup_fitness();

  evaluations_++;
  history_self_.set_fitness(fitness_vec);

  /* OLD
  if (evaluations_ >= 1) {
    // Average fitness values.
    evaluations_++;
    history_self_.avg_fitness(fitness_vec);

  } else {
    evaluations_ = 1;
    history_self_.set_fitness(fitness_vec);
  }
  */


  if (history_self_.get_how_created() == CREATE_MUTATE) {
    gNum_Mutations++;
    result = better(maximize);
    if (result == 1) {
      gNum_Mutations_Better++;
    } else if (result == 0) {
      gNum_Mutations_Equal++;
    }

  } else if (history_self_.get_how_created() == CREATE_RECOMBINE) {
    gNum_Recombinations++;
    result = better(maximize);
    if (result == 1) {
      gNum_Recombinations_Better++;
    } else if (result == 0) {
      gNum_Recombinations_Equal++;
    }
  }

  //  CmnClass::Logger::Msg() << "Op count size: " << history_self_.get_op_count()->size() << endl;
  op_update(*history_self_.get_op_count());
}

void Individual :: set_fitness(vector<double>& fitness_vec)
{
  set_fitness(true, fitness_vec);
}


void Individual :: set_fitness(double fitness)
{
  vector<double> v(1, fitness);
  set_fitness(v);
}


void Individual :: reset_num_evaluations()
{
  evaluations_ = 0;
}


int Individual :: get_num_evaluations()
{
  return evaluations_;
}

void Individual :: incr_num_evaluations()
{
  evaluations_++;
}


int Individual :: get_age() const
{
  return history_self_.get_age();
}

int Individual :: get_age(int evals, int num_individs) const
{
  return history_self_.get_age(evals, num_individs);
}

void Individual :: set_age(int age)
{
  history_self_.set_age(age);
}

void Individual :: increase_age()
{
  history_self_.incr_age();
}



int Individual :: get_age_move() const
{
  return age_move_;
}

void Individual :: set_age_move(int age)
{
  age_move_ = age;
}



void Individual :: set_take_age_older()
{
  assign_age_type_ = Age_Take_Older;
}
void Individual :: set_take_age_average()
{
  assign_age_type_ = Age_Take_Average;
}
void Individual :: set_take_age_younger()
{
  assign_age_type_ = Age_Take_Younger;
}




void Individual :: set_creation(int type)
{
  history_self_.set_how_created(type);
}

int Individual :: get_creation()
{
  return history_self_.get_how_created();
}


void Individual :: print_history()
{
  history_self_.print(cout);
}


unsigned int Individual :: num_features() const
{
  return features_.size();
}

vector<double> Individual :: get_features() const
{
  //  CmnClass::Logger::Msg() << "individual :: get_features() - error, should be instantiated in a subclass." << endl;
  //  while (1) ;
  return features_;
}


void Individual :: make_random()
{
  design_id_ = 0;

  evaluations_ = 0;
  history_self_.clear();
  history_self_.set_how_created(CREATE_RANDOM);

  history_parent1_.clear();
  history_parent2_.clear();
}

void Individual :: duplicate_settings(const Individual *individ2)
{
  assign_age_type_ = individ2->assign_age_type_;
}


void Individual :: duplicate(const Individual *ind_src)
{
  design_id_ = ind_src->design_id_;

  evaluations_ = ind_src->evaluations_;
  history_self_ = ind_src->history_self_;
  assign_age_type_ = ind_src->assign_age_type_;

  history_parent1_ = ind_src->history_parent1_;
  history_parent2_ = ind_src->history_parent2_;

  features_ = ind_src->features_;
}



bool Individual :: mutate(double scale)
{
  //  cout << "Individual :: mutate(d)" << endl;
  design_id_ = 0; // No longer valid.
  evaluations_ = 0;

  history_parent1_ = history_self_;
  history_self_.incr_num_mutate();
  history_self_.set_how_created(CREATE_MUTATE);

  return true;
}

bool Individual :: mutate(double scale, vector<Individual*>&)
{
  //  cout << "Individual :: mutate(d, vector<Ind*>&)" << endl;
  return Individual::mutate(scale);
}

bool Individual :: mutate()
{
  //  cout << "Individual :: mutate()" << endl;
  return Individual::mutate(1.0);
}

bool Individual :: mutate(vector<Individual*>&)
{
  //  cout << "Individual :: mutate(vector<Ind*>&)" << endl;
  return Individual::mutate(1.0);
}




bool Individual :: recombine(double scale, Individual *parent2)
{
  design_id_ = 0; // No longer valid.
  evaluations_ = 0;

  // History for this individual.
  history_parent1_ = history_self_;
  history_parent2_ = parent2->history_self_;
  history_self_.incr_num_recombine();
  history_self_.set_how_created(CREATE_RECOMBINE);  

  return true;
}

bool Individual :: recombine_rand2(double scale, Individual *parent2)
{
  return Individual::recombine(scale, parent2);
}

bool Individual :: recombine(double scale, vector<Individual*>& ind2s)
{
  if (ind2s.size() == 0) {
    CmnClass::Logger::ErrorDTS() << "individual::recombine(vec<Ind*>) - ind2 size == 0." << endl;
    return false;
  }

  return Individual::recombine(scale, ind2s[0]);
}


bool Individual :: recombine(Individual *parent2)
{
  return Individual::recombine(1, parent2);
}

bool Individual :: recombine_rand2(Individual *parent2)
{
  return Individual::recombine(1.0, parent2);
}

bool Individual :: recombine(vector<Individual*>& ind2s)
{
  return Individual::recombine(1.0, ind2s[0]);
}



/*************************************************/

//int Individual :: phenotype_distance(Individual* ind2)
int Individual :: phenotype_distance(Individual*)
{
  // Should be implented in subclass (if appropriate).
  return 0;
}


//bool Individual :: make_phenotype(void* arg)
bool Individual :: make_phenotype(void*)
{
  // Should be implented in subclass (if appropriate).
  return true;
}

bool Individual :: make_phenotype()
{
  return true;
}



double Individual :: distance_apart(const Individual*) const
{
  // Should be implented in subclass (if appropriate).
  CmnClass::Logger::Msg() << "individual :: distance_apart() - should be implemented by subclass." << endl;
  while (1) ;
  return 0;
}


/*************************************************/




///////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////


void Individual :: write_log(char *string)
{
  int loc = 0;

  for (int i=0; i<gNum_Fitness; i++) {
    sprintf(&string[loc], "%.2f %.2f ", history_self_.get_fitness(i),
	    fitness_hist_diff(i));
    while (string[loc] != 0) {
      ++loc;
    }
  }

  if (history_self_.get_how_created() == CREATE_INIT) {
    sprintf(&string[loc], "3\n");
  } else if (history_self_.get_how_created() == CREATE_MUTATE) {
    sprintf(&string[loc], "1\n");
  } else if (history_self_.get_how_created() == CREATE_RECOMBINE) {
    sprintf(&string[loc], "2\n");
  } else if (history_self_.get_how_created() == CREATE_RANDOM) {
    sprintf(&string[loc], "4\n");
  } else {
    sprintf(&string[loc], "-1\n");
  }
}


ostream& Individual :: write_log(ostream& log_stream)
{
  //  printf("write log\n");

  for (int i=0; i<gNum_Fitness; i++) {
    log_stream << history_self_.get_fitness(i) << " ";
  }

  log_stream << get_age() << " ";

  for (int i=0; i<HIST_NUM_COMPLEXITY; i++) {
    log_stream << history_self_.get_complexity(i) << " ";
  }

  if (history_self_.get_how_created() == CREATE_INIT) {
    log_stream << "3";
  } else if (history_self_.get_how_created() == CREATE_MUTATE) {
    log_stream << "1";
  } else if (history_self_.get_how_created() == CREATE_RECOMBINE) {
    log_stream << "2";
  } else if (history_self_.get_how_created() == CREATE_RANDOM) {
    log_stream << "4";
  } else {
    log_stream << "-1";
  }

  if ((history_self_.get_how_created() == CREATE_MUTATE) ||
      (history_self_.get_how_created() == CREATE_RECOMBINE)) {
    for (int i=0; i<gNum_Fitness; i++) {
      log_stream << " " << history_parent1_.get_fitness(i);
    }

    for (int i=0; i<HIST_NUM_COMPLEXITY; i++) {
      log_stream << " " << history_parent1_.get_complexity(i);
    }

  } else {
    for (int i=0; i<gNum_Fitness; i++) {
      log_stream << " 0";
    }

    for (int i=0; i<HIST_NUM_COMPLEXITY; i++) {
      log_stream << " 0";
    }
  }

  log_stream << "\n";

  return log_stream;
}


ostream& Individual :: write_features(ostream& ostr) const
{
  CmnClass::write_vec(ostr, features_);

  return ostr;
}




istream& Individual :: read(istream& istr)
{
  string line;
  vector<string> words;

  // Line 1: Read type and version number:
  getline(istr, line);
  split(line, ' ', words);
  if (words.size() != 2) {
    CmnClass::Logger::Msg() << "individual :: read() - error reading version:" << line << endl;
    while (1) ;
    throw runtime_error("individual :: read() - error reading version:" + line);
  }

  if (words[0] != "IndV") {
    CmnClass::Logger::Msg() << "individual :: read() - error reading version[0]:" << line << endl;
    while (1) ;
    throw runtime_error("individual :: read() - error, expecting IndV, got:" + line);
  }

  int version = from_string<int>(words[1]);
  if (version != INDIVIDUAL_VERSION) {
    CmnClass::Logger::Msg() << "individual :: read() - error reading version[1]:" << line << endl;
    while (1) ;
    throw runtime_error("individual :: read() - error, bad version:" + words[1]);
  }  

  
  // Line 2: Read ID:
  getline(istr, line);
  split(line, ' ', words);
  if (words.size() != 1) {
    CmnClass::Logger::Msg() << "individual :: read() - error parsing: ]" << line << "[" << endl;
    throw runtime_error("individual :: read() - error config line: " + line);
  }
  design_id_.set_str(words[0], 10);


  istr >> history_self_;
  //  istr >> history_parent1_;
  //  istr >> history_parent2_;

  
  return istr;
}


bool Individual :: read(const char *fname)
{
  ifstream file(fname);
  try {
    read(file);
    //  } catch (runtime_error e) {
  } catch (...) {
    return false;
  }

  file.close();
  return true;
}


ostream& Individual :: write(ostream& ostr) const
{
  ostr << "IndV " << INDIVIDUAL_VERSION << endl;
  ostr << design_id_ << endl;
  history_self_.write(ostr);
  //  ostr << history_parent1_;
  //  ostr << history_parent2_;
  return ostr;
}


bool Individual :: write(const char *fname) const
{
  ofstream file(fname);
  write(fname);
  file.close();
  return true;
}



/*************************************************/

void Individual :: print_history_ops()
{
  history_self_.print_ops(cout);
}


void Individual :: print_results()
{
  if (history_self_.get_how_created() == CREATE_MUTATE) {
    if (history_self_.get_fitness() > history_parent1_.get_fitness())
      printf("  Higher!");
    else if (history_self_.get_fitness() == history_parent1_.get_fitness())
      printf("  Same.");
    /*
    else
      printf("  lower. %f < %f", history_self_.get_fitness(),
	     history_parent1_.get_fitness());
    */
 } else if (history_self_.get_how_created() == CREATE_RECOMBINE) {
    if ((history_self_.get_fitness() > history_parent1_.get_fitness()) &&
	(history_self_.get_fitness() > history_parent2_.get_fitness()))
      printf("  Highest!");
    else if ((history_self_.get_fitness() == history_parent1_.get_fitness()) ||
	     (history_self_.get_fitness() == history_parent2_.get_fitness()))
      printf("  Same.");
    else if ((history_self_.get_fitness() > history_parent1_.get_fitness()) ||
	     (history_self_.get_fitness() > history_parent2_.get_fitness()))
      printf("  Middle.");
    /*
    else
      printf("  lower. %f < %f & %f", history_self_.get_fitness(),
	     history_parent1_.get_fitness(), history_parent2_.get_fitness());
    */
  }
}


void Individual :: print_history_full()
{
  streamsize prec = cout.precision();

  if (evaluations_ > 1) {
    double diff, val;

    cout << "Evl(" << evaluations_ << ")                     ";
    val = get_fitness();
    diff = val - history_self_.get_fitness_prev();
    history_parent1_.print_f1(cout);
    cout << " -> ";

    cout << setprecision(6);
    cout << val + diff;
    cout << " diff: " << 2.0*diff;
    //    printf("%6.0f", val + diff);
    //    printf(" diff: %6.0f ", 2.0f*diff);
    return;

  } else if (history_self_.get_how_created() == CREATE_INIT) {
    cout << "Init ";
    //    printf("                     ");
    history_self_.print(cout);
    return;

  } else if (history_self_.get_how_created() == CREATE_RANDOM) {
    cout << "Rand ";
    //    printf("                     ");
    history_self_.print(cout);
    return;

  } else if (history_self_.get_how_created() == CREATE_MUTATE) {
    cout << "Mut ";
    cout << "                         ";
    history_parent1_.print_f1(cout);

  } else if (history_self_.get_how_created() == CREATE_RECOMBINE) {
    cout << "Rec ";
    history_parent1_.print_f1(cout);
    cout << " + ";
    history_parent2_.print_f1(cout);

  }

  cout << " -> ";
  history_self_.print_f1(cout);

  print_results();

  cout.precision(prec);
}


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////

// Returns:
//  true: if ind1 is not as good as ind2
//  false: otherwise
bool individ_compare_max(void const *elt1, void const *elt2)
{
  Individual* ind1 = (Individual*)elt1;
  Individual* ind2 = (Individual*)elt2;

  if (ind1->compare_fitness(true, ind2) != 1) {
    return true;
  }
  return false;
}

bool individ_compare_min(const void *elt1, const void *elt2)
{
  Individual* ind1 = (Individual*)elt1;
  Individual* ind2 = (Individual*)elt2;

  if (ind1->compare_fitness(false, ind2) != 1) {
    return true;
  }
  return false;
}


/***********************************************************/


void individual_print_op_history()
{
  if (gNum_Operators == 0)
    return;

  cout << "m:" << gNum_Mutations_Better << ":" << gNum_Mutations_Equal
       << ":" << gNum_Mutations
       << ", r:" << gNum_Recombinations_Better << ":"
       << gNum_Recombinations_Equal << ":" << gNum_Recombinations;
}


void individual_print_op_history(char *string)
{
  int i, loc;

  if (gNum_Operators == 0)
    return;

  sprintf(string, "[");
  loc = 1;
  for (i=0; i<gNum_Operators-1; i++) {
    sprintf(&string[loc], "%2d ]", gOp_Used_Hist[i]);
    while (string[loc] != ']')
      loc++;
  }
  sprintf(&string[loc], "%2d] ", gOp_Used_Hist[i]);

  while (string[loc] != ' ')
    loc++;

  sprintf(&string[loc], "  m:%d:%d:%d r:%d:%d:%d\n",
	  gNum_Mutations_Better, gNum_Mutations_Equal, gNum_Mutations,
	  gNum_Recombinations_Better, gNum_Recombinations_Equal, gNum_Recombinations);
}


void individual_write_op_history(char *string)
{
  int i, loc;

  sprintf(string, "%d\n", gNum_Operators);
  if (gNum_Operators == 0)
    return;

  loc = 1;
  for (i=0; i<gNum_Operators; i++) {
    sprintf(&string[loc], " %d\n", gOp_Used_Hist[i]);
    while (string[loc] != '\n')
      loc++;
  }

  sprintf(&string[loc], "  m:%d:%d:%d r:%d:%d:%d\n",
	  gNum_Mutations_Better, gNum_Mutations_Equal,
	  gNum_Mutations, gNum_Recombinations_Better,
	  gNum_Recombinations_Equal, gNum_Recombinations);
}


void individual_read_op_history(const char *string)
{
  int i, loc, num;


  num = read_integer(string);
  if (num != gNum_Operators) {
    printf("individual_read_op_history :: # operators does not match: %d != %d.\n",
	   num, gNum_Operators);
  }

  loc = 0;
  for (i=0; i<num; i++) {
    while (string[loc] != ' ') {
      loc++;
      if (loc > 1000) {
	printf("individual_read_op_history :: error reading op hist %d:%s\n",
	       loc, string);
	return;
      }
    }
    loc++;
    gOp_Used_Hist[i] = read_integer(&string[loc]);
  }

  while (string[loc] != 'm') {
    loc++;
    if (loc > 1000) {
      printf("individual_read_op_history :: error finding varation history %d:%s\n",
	     loc, string);
      return;
    }
  }

  num = sscanf(&string[loc], "m:%d:%d:%d r:%d:%d:%d\n",
	       &gNum_Mutations_Better, &gNum_Mutations_Equal,
	       &gNum_Mutations, &gNum_Recombinations_Better,
	       &gNum_Recombinations_Equal, &gNum_Recombinations);
  if (num != 6) {
    printf("individual_read_op_history :: error read %d:  m:%d:%d:%d r:%d:%d:%d\n",
	   num, gNum_Mutations_Better, gNum_Mutations_Equal,
	   gNum_Mutations, gNum_Recombinations_Better,
	   gNum_Recombinations_Equal, gNum_Recombinations);
  }
}
