/***

This file is part of the OptAlgs C++ library originally developed
by Gregory Hornby.

This code is released under the BSD-3 license:
  http://opensource.org/licenses/BSD-3-Clause
See the file "license.txt" in the root directory for full details.

***/


/*

 author: Gregory S. Hornby
 file: opt_alg.h

*/


#ifndef OPTALG_HEADER_FILE
#define OPTALG_HEADER_FILE

#include <iostream>
#include <fstream>
#include <map>
#include <semaphore.h>
#include <set>
#include <string>
#include <vector>


class FitnessFunc;
class Individual;

class OptAlg
{
  friend std::ostream& operator<<(std::ostream&, const OptAlg&);
  friend std::istream& operator>>(std::istream&, OptAlg&);

public:
  static bool is_print_log_;
  const static std::string class_name_;

  static std::vector<Individual*> individual_types_;
  static void add_individual_type(const Individual*);
  static const Individual* get_individual_type(const std::string& type_name);
  static void clear_individual_types();

  static std::vector<OptAlg*> alg_types_;
  static void add_alg_type(const OptAlg*);
  static const OptAlg* get_alg_type(const std::string& type_name);
  static void clear_alg_types();

  static void clear_types();
  static void set_random_seed(unsigned long seed);
  

  /* ******* */
  OptAlg(const FitnessFunc* func, const Individual* ind_template);
  virtual ~OptAlg();

  virtual std::string get_class_name() const { return class_name_; }
  const Individual* get_ind_template() const { return ind_template_; }

  virtual void copy_settings(const OptAlg* p_src);

  virtual OptAlg* new_instance() const = 0;
  virtual OptAlg* make_copy() const = 0;

  virtual void clear();
  virtual void init_variables(); // Only resets evaluation counts.


  void make_filename(const char *ending, char *fname) const;
  int set_logfile(const char *n);

  void set_print_debug(bool f);
  void set_print_gen_stats(bool f);
  void set_save_best(bool b);


  bool get_maximize() const { return is_maximizing_; };
  void set_maximize() { is_maximizing_ = true; };
  void set_minimize() { is_maximizing_ = false; };

  virtual int configure(const char *fname, int verbose);


  // **************************************************

  unsigned int get_num_evals() const;
  void increment_num_evals(unsigned int num);
  void increment_num_evals();

  virtual void do_search_step() = 0;


  // **************************************************

  virtual void print_current_individs(std::ostream& ostr) = 0;
  virtual void print_current_individs() = 0;


  bool read_sample_individ(const char *fname);
  bool read_seed_individ(const char *fname);

  virtual std::ostream& write_template(std::ostream& ostr) const;
  virtual std::ostream& write_header(std::ostream& ostr) const;
  virtual std::ostream& write(std::ostream& os) const;
  virtual bool write(const char *fname) const;
  virtual bool write_backup() const;

  virtual std::istream& read(std::istream& is);
  virtual bool read(const char *fname);


  // *********************


protected:
  bool is_maximizing_;
  bool is_print_debug_info_;
  bool is_run_sample_;
  bool is_save_best_;
  bool is_save_log_;

  // Generated statistics:
  unsigned int num_evaluations_;
  double best_fitness_;

  std::string log_file_name_;
  std::ofstream log_file_;

  Individual* ind_template_; // Template individual for configuring things, like number of genes/inputs.
  Individual* ind_seed_; // Seed for an evolutionary run (supercedes Sample_Individ).

  FitnessFunc* fitness_func_;
};

#endif
