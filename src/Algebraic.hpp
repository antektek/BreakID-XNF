#pragma once

#include "global.hpp"

class Specification;
class Breaker;

class Permutation : public std::enable_shared_from_this<Permutation> {
private:
  std::unordered_map<unsigned int, unsigned int> perm;
  std::vector<unsigned int> cycleReprs; // smallest lit in each cycle
  unsigned int maxCycleSize;
  size_t hash;

  

public:
  std::vector<unsigned int> domain;
  std::vector<unsigned int> posDomain;
  std::vector<unsigned int> image;
  
  void addFromTo(unsigned int from, unsigned int to);
  void addCycle(std::vector<unsigned int>& cyc);

  Permutation();
  Permutation(std::vector<std::pair<unsigned int, unsigned int> >& tuples);
  // Permutation constructed from swapping two rows.
  Permutation(std::vector<unsigned int>& row1, std::vector<unsigned int>& row2);

  ~Permutation() {
  };

  unsigned int getImage(unsigned int from);
  // return value is true iff the image is different from the original
  bool getImage(std::vector<unsigned int>& orig, std::vector<unsigned int>& img);
  bool getImage(std::map<unsigned int,unsigned int>& orig, std::map<unsigned int,unsigned int>& img);
  void getCycle(unsigned int lit, std::vector<unsigned int>& orb);
  bool isInvolution();
  bool permutes(unsigned int lit);
  unsigned int supportSize();
  bool isIdentity();

  void print(std::ostream& out);
  void printAsProjection(std::ostream& out);

  bool formsMatrixWith(sptr<Permutation> other);
  std::pair<sptr<Permutation>, sptr<Permutation> > getLargest(sptr<Permutation> other);
  void getSharedLiterals(sptr<Permutation> other, std::vector<unsigned int>& shared);
  std::vector<unsigned int>& getCycleReprs();
  unsigned int getMaxCycleSize();
  unsigned int getNbCycles();

  bool equals(sptr<Permutation> other);
};

class Matrix {
private:
  std::vector<std::vector<unsigned int>* > rows; // TODO: refactor this as 1 continuous vector
  std::unordered_map<unsigned int, unsigned int> rowco;
  std::unordered_map<unsigned int, unsigned int> colco;

public:
  Matrix();
  ~Matrix();
  void print(std::ostream& out);

  void add(std::vector<unsigned int>* row);
  unsigned int nbColumns();
  unsigned int nbRows();
  void tryToAddNewRow(sptr<Permutation> p, unsigned int rowIndex, Specification* theory);
  std::vector<unsigned int>* getRow(unsigned int rowindex);
  bool permutes(unsigned int x);
  unsigned int getLit(unsigned int row, unsigned int column);

  unsigned int getRowNb(unsigned int x);
  unsigned int getColumnNb(unsigned int x);

  sptr<Permutation> testMembership(const sptr<Permutation> p);
  sptr<Permutation> getProductWithRowsWap(const sptr<Permutation> p, unsigned int r1, unsigned int r2); // return p*swap(r1,r2)
};

class Group {
private:
  std::vector<sptr<Permutation> > permutations;
  std::vector<sptr<Matrix> > matrices;
  std::unordered_set<unsigned int> support;

  void cleanPermutations(sptr<Matrix> matrix); // remove permutations implied by the matrix

public:
  // NOTE: if a group has a shared pointer to a theory, and a theory a shared pointer to a group, none of the memory pointed to by these pointers will ever be freed :(
  Specification* theory; // non-owning pointer

  Group() {
  };

  ~Group() {
  };

  void add(sptr<Permutation> p);
  void checkColumnInterchangeability(sptr<Matrix> m);

  void print(std::ostream& out);

  sptr<Matrix> getInitialMatrix();
  
  void addMatrices();
  void addMatrix(sptr<Matrix> m); // cnf-parameter, otherwise we have to store a pointer to the cnf here :(
  unsigned int getNbMatrices();
  unsigned int getNbRowSwaps();
  
  sptr<Matrix> getMatrix(unsigned int idx);

  void getDisjointGenerators(std::vector<sptr<Group> >& subgroups);
  unsigned int getSize();

  bool permutes(unsigned int lit);
  unsigned int getSupportSize();

  void getOrderAndAddBinaryClausesTo(Breaker& brkr, std::vector<unsigned int>& out_order); // returns a vector containing a lit for literals relevant to construct sym breaking clauses
  void addBinaryClausesTo(Breaker& brkr, std::vector<unsigned int>& out_order, const std::unordered_set<unsigned int>& excludedLits);
  void addBreakingClausesTo(Breaker& brkr);

  void maximallyExtend(sptr<Matrix> matrix, unsigned int indexOfFirstNewRow);
};
