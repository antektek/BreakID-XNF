/* Graph is a class used as a wrapper for saucy methods
 */

#pragma once

#include "global.hpp"
#include "Algebraic.hpp"

// TODO: make a die function for mem-assignment issues?
/*static void die(string& s){
  std::clog << s << endl;
  exit(1);
}*/

class Matrix;
class Group;
struct saucy_graph;

class Graph : public std::enable_shared_from_this<Graph> {
private:

public:
  saucy_graph* saucy_g;
  std::vector<unsigned int> colorcount; // keeps track of the number of times a color is used, so that no color is never used (seems to give Saucy trouble)
  // @INVAR: for all x: colorcount[x]>0

  Graph(std::unordered_set<sptr<Clause>, UVecHash, UvecEqual>& clauses);
  Graph(std::unordered_set<sptr<Clause>, UVecHash, UvecEqual>& CNFClauses, std::unordered_set<sptr<Clause>, UVecHash, UvecEqual>& XORClauses);
  Graph(std::unordered_set<sptr<Rule>, UVecHash, UvecEqual>& rules);
  Graph(std::unordered_set<sptr<PBConstraint>, UVecHash, UvecEqual>& constrs);
  ~Graph();

private:
    //Interaction with saucy:
    void initializeGraph(unsigned int nbNodes, unsigned int nbEdges, std::map<unsigned int, unsigned int> &lit2color, std::vector<std::vector<unsigned int> > &neighbours);
    void freeGraph();
    unsigned int getNbNodesFromGraph();
    unsigned int getNbEdgesFromGraph();
    unsigned int getColorOf(unsigned int node);
    unsigned int nbNeighbours(unsigned int node);
    unsigned int getNeighbour(unsigned int node, unsigned int nbthNeighbour);
    void setNodeToNewColor(unsigned int node);
    void getSymmetryGeneratorsInternal(std::vector<sptr<Permutation> > &out_perms);


public:
  unsigned int getNbNodes();
  unsigned int getNbEdges();
  void print();
  void setUniqueColor(unsigned int lit);
  void setUniqueColor(const std::vector<unsigned int>& lits);
  void getSymmetryGenerators(std::vector<sptr<Permutation> >& out_perms);
};

