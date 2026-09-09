#pragma once

#include "global.hpp"

class Permutation;
class Specification;

class Breaker {
private:
  std::unordered_set<sptr<Clause>, UVecHash, UvecEqual> clauses;
  std::unordered_set<sptr<PBConstraint>, UVecHash, UvecEqual> pbConstrs;
  sptr<Specification> originalTheory;
  unsigned int nbExtraVars = 0;
  unsigned int nbBinClauses = 0;
  unsigned int nbRowClauses = 0;
  unsigned int nbRegClauses = 0;

  void addUnary(unsigned int l1);
  void addBinary(unsigned int l1, unsigned int l2);
  void addTernary(unsigned int l1, unsigned int l2, unsigned int l3);
  void addQuaternary(unsigned int l1, unsigned int l2, unsigned int l3, unsigned int l4);
  void add(sptr<Clause> cl);
  void add(sptr<Permutation> perm, std::vector<unsigned int>& order, bool limitExtraConstrs);
  void addPB(sptr<Permutation> perm, std::vector<unsigned int>& order, bool limitExtraConstrs);
  void addShatter(sptr<Permutation> perm, std::vector<unsigned int>& order, bool limitExtraConstrs);

public:
  Breaker(sptr<Specification> origTheo);

  ~Breaker() {
  };

  //Prints the current breaker
  void print();

  void addBinClause(unsigned int l1, unsigned int l2);
  void addRegSym(sptr<Permutation> perm, std::vector<unsigned int>& order);
  void addRowSym(sptr<Permutation> perm, std::vector<unsigned int>& order);

  unsigned int getAuxiliaryNbVars();
  unsigned int getTotalNbVars();
  unsigned int getAddedNbClauses();
  unsigned int getTotalNbClauses();

  unsigned int getNbBinClauses();
  unsigned int getNbRowClauses();
  unsigned int getNbRegClauses();

  unsigned int getTseitinVar();
};
