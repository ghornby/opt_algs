/***

This file is part of the OptAlgs C++ library originally developed
by Gregory Hornby.

This code is released under the BSD-3 license:
  http://opensource.org/licenses/BSD-3-Clause
See the file "license.txt" in the root directory for full details.

***/


/*******************************************************

  File: oaf_rosen.cpp

********************************************************/

#include "oa_frosen.h"

#include "individ_real.h"


const std::string OAFRosen :: class_name_("OAF_Rosen");


OAFRosen :: OAFRosen()
{
  num_evaluations_ = 0;
}

OAFRosen :: ~OAFRosen()
{
}


FitnessFunc* OAFRosen :: new_instance() const
{
  FitnessFunc* func = new OAFRosen();
  return func;
}


double OAFRosen :: evaluate(Individual* ind)
{
  Individ_Real* ind_real = (Individ_Real*)ind;


  std::vector<double> x(ind_real->get_genes());

  std::vector<double> y = ind_real->get_genes();


  double fitness = 0;
  for (std::vector<double>::size_type i=0; i<x.size(); i++) {
    double val = 1.0 - x[i];
    fitness += val * val;
  }

  for (std::vector<double>::size_type i=1; i<x.size(); i++) {
    double val = x[i] - x[i-1]*x[i-1];
    fitness += 100.0 * val * val;
  }

  ind_real->set_fitness(fitness);
  update_best(fitness);
  increment_num_evals();

  return fitness;
}


