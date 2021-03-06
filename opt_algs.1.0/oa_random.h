/***

This file is part of the OptAlgs C++ library originally developed
by Gregory Hornby.

This code is released under the BSD-3 license:
  http://opensource.org/licenses/BSD-3-Clause
See the file "license.txt" in the root directory for full details.

***/


/*

 author: Gregory S. Hornby
 file: oa_random.h

 Description:
 This a simple random search algorithm.

*/


#ifndef OA_RANDOM_HEADER_FILE
#define OA_RANDOM_HEADER_FILE

#include "opt_alg.h"

#include <iostream>
#include <fstream>
#include <map>
#include <semaphore.h>
#include <set>
#include <string>
#include <vector>


class FitnessFunc;
class Individual;

class OA_Random : public OptAlg
{
  friend std::ostream& operator<<(std::ostream&, const OA_Random&);
  friend std::istream& operator>>(std::istream&, OA_Random&);

public:
  const static std::string class_name_;


  /* ******* */
  OA_Random(const FitnessFunc* func, const Individual* ind_template);
  ~OA_Random();

  std::string get_class_name() const { return class_name_; }

  void copy_settings(const OA_Random* p_src);

  OptAlg* new_instance(const Individual*) const;
  OptAlg* new_instance() const;
  OptAlg* make_copy() const;

  void clear();
  void init_variables(); // Only resets evaluation counts.

  virtual int configure(const char *fname, int verbose);


  // **************************************************


  void do_search_step();


  // **************************************************

  void print_current_individs(std::ostream& ostr);
  void print_current_individs();


  virtual std::ostream& write_template(std::ostream& ostr) const;
  virtual std::ostream& write_header(std::ostream& ostr) const;
  virtual std::ostream& write(std::ostream& os) const;
  virtual bool write(const char *fname) const;
  virtual bool write_backup() const;

  virtual std::istream& read(std::istream& is);
  virtual bool read(const char *fname);


  // *********************


protected:
  Individual* individ_best_;
  Individual* individ_new_;
};

#endif
