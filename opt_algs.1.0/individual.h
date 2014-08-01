/***

This file is part of the OptAlgs C++ library originally developed
by Gregory Hornby.

This code is released under the BSD-3 license:
  http://opensource.org/licenses/BSD-3-Clause
See the file "license.txt" in the root directory for full details.

***/


#ifndef INDIVIDUAL_HEADER_FILE
#define INDIVIDUAL_HEADER_FILE

#include <iostream>
#include <fstream>
#include <gmpxx.h>

#include "ind_history.h"

const double Worst_Fitness = 1.5E+15;
const int NUM_HISTORY_OPS = 15;
const int NUM_FITNESS = 20;
const int HIST_NUM_COMPLEXITY = 27;


extern int gNum_Fitness;
extern int gNum_Operators;
extern double gRec_Delta;
extern double gMutate_Size;


bool op_incr_used(int op_num);
int op_get_last(void);
void op_reset(void);
void op_update(std::vector<int>& operation);


int recombine_int(int val1, int val2);
double recombine_double_as_int(double val1, double val2);
double recombine_double(double val1, double val2);
double recombine(double val1, double val2);
void recombine_elliptic(std::vector<double>& p1, std::vector<double>& p2,
			std::vector<double>& c1);



class Individual
{
  friend std::ostream& operator<<(std::ostream&, const Individual&);
  friend std::istream& operator>>(std::istream&, const Individual&);

public:
  static bool is_print_log_;
  const static std::string class_name_;

  Individual();
  Individual(const Individual *sample);
  virtual ~Individual();

  bool is_deleted() const;
  virtual std::string get_class_name() const = 0;
  virtual Individual* new_instance() const = 0;


  /*** Interface Methods ***/
  virtual void clear();
  virtual bool valid() const;

  void set_id(const mpz_class& id) { design_id_ = id; };
  mpz_class get_id() const { return design_id_; };
  std::string get_id_string() const { return design_id_.get_str(); };


  double fitness_hist_diff(int index) const;
  double fitness_hist_diff() const;
  virtual int better(bool is_maximizing) const;
  virtual bool is_keep(bool is_maximize) const;
  virtual int compare_fitness(bool is_maximizing, Individual *individ2);
  virtual bool same(const Individual *individ2) const;
  bool created_by_variation() const;

  double get_fitness() const;
  void get_fitness(std::vector<double>& fitness_vec) const;
  const std::vector<double>& get_fitness_vec() const;
  void set_fitness(double fitness);
  void set_fitness(bool maximize, const std::vector<double>& fitness_vec);
  void set_fitness(std::vector<double>& fitness_vec);

  void set_eval_ok() { is_eval_ok_ = true; }
  void set_eval_fail() { is_eval_ok_ = false; }

  void reset_num_evaluations();
  int get_num_evaluations() const;
  void incr_num_evaluations();

  int get_age() const;
  int get_age(int evals, int num_individs) const;
  void set_age(int age);
  void increase_age();

  int get_age_move() const;
  void set_age_move(int age);

  void set_take_age_older();
  void set_take_age_average();
  void set_take_age_younger();

  void set_creation(int type);
  int get_creation() const;
  
  void print_history();

  virtual unsigned int num_features() const;
  virtual std::vector<double> get_features() const;

  virtual void make_random();
  virtual void duplicate_settings(const Individual *individ2);
  virtual void duplicate(const Individual *ind);

  virtual bool mutate(double scale);
  virtual bool mutate(double scale, std::vector<Individual*>& ind2s);
  virtual bool mutate();
  virtual bool mutate(std::vector<Individual*>& ind2s);



  virtual bool recombine(double scale, Individual *parent2);
  virtual bool recombine_rand2(double scale, Individual *parent2);
  virtual bool recombine(double scale, std::vector<Individual*>& ind2s);

  virtual bool recombine(Individual *parent2);
  virtual bool recombine_rand2(Individual *parent2);
  virtual bool recombine(std::vector<Individual*>& ind2s);


  virtual int phenotype_distance(Individual* ind2);
  virtual bool make_phenotype(void* arg);
  virtual bool make_phenotype();
  virtual double distance_apart(const Individual* ind2) const;


  void write_log(char *string);
  std::ostream& write_log(std::ostream& log_stream);

  virtual std::ostream& write_features(std::ostream& ostr) const;

  virtual std::istream& read(std::istream& istr);
  virtual bool read(const char *fname);
  virtual std::ostream& write(std::ostream& ostr) const;
  virtual bool write(const char *fname) const;

  void print_history_ops();
  void print_results();
  void print_history_full();


protected:
  std::vector<double> features_;


private:
  bool is_deleted_;
  bool is_eval_ok_;
  mpz_class design_id_;
  unsigned int num_evaluations_; // Count of number of times evaluated.
  int age_move_; // For steady-state ALPS: try_move_up()
  int assign_age_type_;


  IndHistory history_self_;
  IndHistory history_parent1_;
  IndHistory history_parent2_;
};

bool individ_compare_max(void const *elt1, void const *elt2);
bool individ_compare_min(const void *elt1, const void *elt2);

void individual_print_op_history();
void individual_print_op_history(char *string);
void individual_write_op_history(char *string);
void individual_read_op_history(const char *string);


#endif
