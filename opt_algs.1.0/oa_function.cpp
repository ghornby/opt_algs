/***

This file is part of the OptAlgs C++ library originally developed
by Gregory Hornby.

This code is released under the BSD-3 license:
  http://opensource.org/licenses/BSD-3-Clause
See the file "license.txt" in the root directory for full details.

***/


/*******************************************************

  File: oa_function.cpp

********************************************************/

#include "oa_function.h"

using namespace std;


const std::string FitnessFunc :: class_name_("FitnessFunc");


FitnessFunc :: FitnessFunc()
  : is_minimizing_(true), num_evaluations_(0), best_fitness_(0)
{
}

FitnessFunc :: ~FitnessFunc()
{
}


void FitnessFunc :: update_best(double val)
{
  if (num_evaluations_ == 0) {
    best_fitness_ = val;
    return;
  }

  if (is_minimizing_) {
    if (val < best_fitness_) {
      best_fitness_ = val;
    }

  } else {
    if (val > best_fitness_) {
      best_fitness_ = val;
    }
  }
}
