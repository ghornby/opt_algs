/***

This file is part of the OptAlgs C++ library originally developed
by Gregory Hornby.

This code is released under the BSD-3 license:
  http://opensource.org/licenses/BSD-3-Clause
See the file "license.txt" in the root directory for full details.

***/

/*******************************************************

  File: main.cpp

********************************************************/

#include "individ_real.h"
#include "oa_frosen.h"
#include "oa_hillclimb.h"
#include "oa_random.h"


#include "cmn_random.h"
#include "cmn_string_proc.h"
#include "cmn_utils.h"



#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

#include <sys/time.h>


/***********************************/


int main(int argc, char **argv)
{
  cout << "OA Example 1." << endl;

  CmnClass::seed_random(0);


  const unsigned int Number_Genes = 20;
  Individ_Real* individ_config = new Individ_Real(Number_Genes);

  FitnessFunc* func = new OAFRosen();
  
  //  OptAlg* oa = new OA_Random(func, individ_config);
  OptAlg* oa = new HillClimber(func, individ_config);
  
  for (unsigned int i=0; i<200000; i++) {
    oa->do_search_step();
  }

  return 0;
}


