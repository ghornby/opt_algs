/***

This file is part of the OptAlgs C++ library originally developed
by Gregory Hornby.

This code is released under the BSD-3 license:
  http://opensource.org/licenses/BSD-3-Clause
See the file "license.txt" in the root directory for full details.

***/


/************************************************************
 author: Gregory S. Hornby.
************************************************************/
#ifndef OA_FUNCTION_HEADER_FILE
#define OA_FUNCTION_HEADER_FILE

#include <iostream>
#include <fstream>
#include <vector>

class Individual;

class FitnessFunc
{
 public:
  FitnessFunc();
  virtual ~FitnessFunc();

  virtual std::string get_class_name() const { return class_name_; }

  virtual FitnessFunc* new_instance() const = 0;
  virtual double evaluate(Individual*) = 0;

  void increment_num_evals();
  unsigned int get_num_evaluations() const { return num_evaluations_; };
  double get_best_fitness() const { return best_fitness_; };

 protected:
  void update_best(double);

  bool is_minimizing_;
  
  const static std::string class_name_;

  unsigned int num_evaluations_;
  double best_fitness_;
};

#endif
