/***

This file is part of the OptAlgs C++ library originally developed
by Gregory Hornby.

This code is released under the BSD-3 license:
  http://opensource.org/licenses/BSD-3-Clause
See the file "license.txt" in the root directory for full details.

***/


/************************************************************
  Copyright Gregory S. Hornby.

************************************************************/
#ifndef INDIVID_REAL_HEADER_FILE
#define INDIVID_REAL_HEADER_FILE

#include "individual.h"
#include <vector>


class Individ_Real : public Individual
{
  friend std::ostream& operator<<(std::ostream&, const Individ_Real&);
  friend std::istream& operator>>(std::istream&, Individ_Real&);


public:
  const static std::string class_name_;

  /********** Methods **********/
  Individ_Real(int num_genes = 0); // perhaps input data struct
  Individ_Real(const Individual *sample);
  ~Individ_Real();

  /*** Interface Methods ***/
  std::string get_class_name() const { return class_name_; }
  Individ_Real* new_instance() const;
  void clear();
  bool same(const Individ_Real *individ2);

  void set_minmax(std::vector<double>& min,
		  std::vector<double>& max);
  void set_init_minmax(std::vector<double>& min,
		       std::vector<double>& max);
  void set_minmax(double min, double max);
  void set_init_minmax(double min, double max);
  void zero_genes();
  int get_num_genes() const;
  bool set_num_genes(int num_genes);
  const std::vector<double>& get_genes() const;

  void set_genes(std::vector<double>& genes);

  unsigned int num_features() const;
  std::vector<double> get_features() const;
  std::ostream& write_features(std::ostream& ostr) const;


  void make_random();
  void duplicate_settings(const Individ_Real *individ2);
  void duplicate_settings(const Individual *individ2);
  void duplicate(const Individ_Real *ind);
  void duplicate(const Individual *ind);

  bool mutate(double scale);
  bool mutate();
  void get_unique(std::vector<Individual*>& inds2,
		  std::vector<Individ_Real*>& indreals);
  bool mutate(double scale, std::vector<Individual*>& inds2);
  bool mutate(std::vector<Individual*>& inds2);

  void recombine_linear_clip(double dist, std::vector<double>& parent1,
			     std::vector<double>& parent2);
  void recombine_linear(double dist, std::vector<double>& parent1,
			std::vector<double>& parent2);

  bool recombine_rand2(double scale, Individ_Real *parent2);
  bool recombine_rand2(double scale, Individual *parent2);
  bool recombine_rand2(Individ_Real *parent2);
  bool recombine_rand2(Individual *parent2);

  bool recombine(double scale, Individ_Real *parent2);
  bool recombine(Individ_Real *parent2);

  bool recombine(double scale, Individual *parent2);
  bool recombine(Individual *parent2);

  bool recombine(double scale, std::vector<Individual*>& inds2);
  bool recombine(std::vector<Individual*>& inds2);

  bool make_phenotype(void* arg);
  bool make_phenotype();
  int phenotype_distance(Individ_Real* ind2);
  int phenotype_distance(Individual* ind2);

  double distance_apart(const Individ_Real* ind2) const;
  double distance_apart(const Individual* ind2) const;


  std::istream& read(std::istream& file);
  bool read(const char *fname);
  std::ostream& write(std::ostream& file) const;
  bool write(const char *fname) const;


 protected:
  int num_genes_;
  std::vector<double> gene_vec_;
  std::vector<double> max_val_;
  std::vector<double> min_val_;
  std::vector<double> init_max_;
  std::vector<double> init_min_;
};


#endif
