#include "global.hpp"
#include "limits.h"
#include <fstream>
#include <stdlib.h>



using namespace std;

unsigned int nVars = 0;
std::vector<unsigned int> fixedLits;
std::string inputFile = "";
std::string inputSymFile = "";
std::string generatorFile = "";
std::string logfile = "";
std::ofstream logstream;

time_t startTime;


// OPTIONS:
bool logging = false;
bool useMatrixDetection = true;
bool useBinaryClauses = true;
bool onlyPrintBreakers = false;
bool useShatterTranslation = false;
bool useFullTranslation = false;
int symBreakingFormLength = 50;
bool aspinput = false;
int nbPBConstraintsToGroup = -1;
bool xnfinput = false;
unsigned int verbosity = 1;
int timeLim = INT_MAX;
unsigned int nbConstraintsInProofLog = 0;

int timeLeft() {
  time_t now;
  time(&now);
  return timeLim - difftime(now, startTime);
}

bool timeLimitPassed() {
  return timeLeft() <= 0;
}

void gracefulError(string str) {
  std::cerr << str << "\nExiting..." << endl;
  exit(1);
}

//Prints a clause c a rule falsevar <- (not c)
//afterwards, falsevar (a new variable) will be added to the "false" constraints.
void Clause::printAsRule(std::ostream& ostr, unsigned int falsevar) {
  ostr << "1 " << falsevar << " ";
  std::set<unsigned int> posBodyLits;
  std::set<unsigned int> negBodyLits;
  for (auto lit : lits) {
    auto decoded = decode(lit);
    if (decoded > 0) {
      posBodyLits.insert(decoded);
    } else {
      negBodyLits.insert(-decoded);
    }
  }
  ostr << (posBodyLits.size() + negBodyLits.size()) << " " << negBodyLits.size() << " ";
  for (auto decodedlit : negBodyLits) {
    ostr << decodedlit << " ";
  }
  for (auto decodedlit : posBodyLits) {
    ostr << decodedlit << " ";
  }
  ostr << std::endl;
}