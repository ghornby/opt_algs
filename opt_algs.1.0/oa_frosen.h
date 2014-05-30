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
#ifndef OAF_ROSENBROCK_HEADER_FILE
#define OAF_ROSENBROCK_HEADER_FILE

#include "oa_function.h"

#include <iostream>
#include <fstream>
#include <vector>


class Individual;

class OAFRosen : public FitnessFunc
{
 public:
  OAFRosen();
  ~OAFRosen();

  std::string get_class_name() const { return class_name_; }

  FitnessFunc* new_instance() const;
  double evaluate(Individual*);

  unsigned int get_num_evaluations() const { return num_evaluations_; };


 protected:
  const static std::string class_name_;
};

#endif
