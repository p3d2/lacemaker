/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(soft_exclude,PairSoftExclude);
// clang-format on
#else

#ifndef LMP_PAIR_SOFT_EXCLUDE_H
#define LMP_PAIR_SOFT_EXCLUDE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSoftExclude : public Pair {
  friend class Pair;

 public:
  PairSoftExclude(class LAMMPS *);
  ~PairSoftExclude();

  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);

 protected:
  double cut_global;
  double **prefactor;
  double **cut;
  int N_ignore;  // new member variable declaration

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif