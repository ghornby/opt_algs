/***

This file is part of the OptAlgs C++ library originally developed
by Gregory Hornby.

This code is released under the BSD-3 license:
  http://opensource.org/licenses/BSD-3-Clause
See the file "license.txt" in the root directory for full details.

***/

#ifndef IND_HISTORY_HEADER_FILE
#define IND_HISTORY_HEADER_FILE

#include <iostream>
#include <fstream>
#include <vector>


const int CREATE_NOT = -1;

const int CREATE_MUTATE = 1;
const int CREATE_RECOMBINE = 2;
const int CREATE_INIT = 3;

const int CREATE_WAITING = 4;
const int CREATE_RANDOM = 5;
const int CREATE_MUTATION = 6;
const int CREATE_RECOMBINATION = 7;
const int CREATE_DUPLICATE = 8;



class IndHistory
{
  friend std::ostream& operator<<(std::ostream&, const IndHistory&);
  friend std::istream& operator>>(std::istream&, IndHistory&);

 public:
  IndHistory();
  ~IndHistory();

  void clear();


  void backup_fitness();
  void avg_fitness(const std::vector<double>& fvec);
  void set_fitness(const std::vector<double>& fvec);

  double get_fitness() const;
  double get_fitness(int index) const;
  const std::vector<double>& get_fitness_vec() const;

  double get_fitness_prev() const;

  double get_complexity(int index) const;


  int get_how_created() const;
  void set_how_created(int h);

  int get_age() const;
  int get_age(int evals, int denom) const;
  void set_age(int a);
  void incr_age();
  void take_older(const IndHistory& h2);
  void take_younger(const IndHistory& h2);
  void take_average(const IndHistory& h2);

  void incr_num_mutate();
  void incr_num_recombine();

  std::vector<int>* get_op_count();


  void print_created_orig(std::ostream& ofs);
  void print_created(std::ostream& ofs);
  void print_ops(std::ostream& ofs);
  void print_fitnesses(std::ostream& ofs);
  void print(std::ostream& ofs);
  void print_f1(std::ostream& ofs);
  std::istream& read(std::istream& istr);
  std::ostream& write(std::ostream& ostr) const;

 private:
  int age_;
  int how_created_;
  int last_op_;
  int num_mutate_;
  int num_recombine_;

  std::vector<int> op_count_;
  std::vector<double> complexity_;
  std::vector<double> fitness_;
  std::vector<double> fitness_prev_;
};

#endif
