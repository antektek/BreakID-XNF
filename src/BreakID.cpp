#include <fstream>
#include <sstream>
#include <iterator>
#include <cstring>
#include <stdlib.h>

#include "Graph.hpp"
#include "global.hpp"
#include "Algebraic.hpp"
#include "Theory.hpp"
#include "Breaking.hpp"

using namespace std;

void setFixedLits(std::string& filename) {
  if (verbosity > 0) {
    std::clog << "*** Reading fixed variables: " << filename << std::endl;
  }

  ifstream file(filename);
  if (!file) {
    gracefulError("No fixed variables file found.");
  }
  
  fixedLits.clear();
  string line;
  while (getline(file, line)) {
    istringstream iss(line);
    int l;
    while (iss >> l) {
      fixedLits.push_back(encode(l));
      fixedLits.push_back(encode(-l));
    }
  }
}

void addInputSym(sptr<Group> group) {
  if(inputSymFile == ""){
  	return; // no input symmetry
  }
  if (verbosity > 0) {
    std::clog << "*** Reading input symmetry from: " << inputSymFile << std::endl;
  }
  
  ifstream file(inputSymFile);
  if (!file) {
    gracefulError("No input symmetry file found.");
  }
  
  string line;
  std::vector<unsigned int> cycle;
  while (getline(file, line)) {
    if (line.front() == '(') { // this is a symmetry generator line    
      if(verbosity>1){
        std::clog << "**** adding input symmetry generator" << std::endl;
      }

      sptr<Permutation> perm = std::make_shared<Permutation>();
      
      stringstream sym(line);
      string cycle_str;

      getline(sym, cycle_str, '('); // skip first empty line
      while(std::getline(sym, cycle_str, '(')){
        stringstream cycle_ss(cycle_str);
        cycle.clear();
        int l=0;
        while(cycle_ss >> l){
          cycle.push_back(encode(l));
        }
        perm->addCycle(cycle);
      }
      group->add(perm);
      perm->print(std::cerr);
    }else if (line.front() == 'r'){ // this is an interchangeability matrix block
      if(verbosity>1){
        std::clog << "**** adding input symmetry matrix" << std::endl;
      }
      unsigned int nbRows;
      unsigned int nbColumns;
      char tmp;
      istringstream iss(line);
      for(unsigned int i=0; i<4; ++i){ iss >> tmp; } // skip "rows"
      iss >> nbRows;
      for(unsigned int i=0; i<7; ++i){ iss >> tmp; } // skip "columns"
      iss >> nbColumns;
      sptr<Matrix> mat = std::make_shared<Matrix>();
      
      for(unsigned int r=0; r<nbRows; ++r){
        std::vector<unsigned int>* newRow = new std::vector<unsigned int>(nbColumns);
        getline(file, line);
        istringstream iss_row(line);
        int tmp;
        for(unsigned int c=0; c<nbColumns; ++c){
          iss_row >> tmp;
          (*newRow)[c]=encode(tmp);
        }
        mat->add(newRow);
      }
      group->addMatrix(mat);
    }else{
      gracefulError("Unexpected line in input symmetry file: "+line);
    }
  }
}

namespace options {
  // option strings:
  string instancefile = "-f";
  string nointch = "-no-row";
  string nobinary = "-no-bin";
  string formlength = "-s";
  string verbosity = "-v";
  string timelim = "-t";
  string help = "-h";
  string nosmall = "-no-small";
  string norelaxed = "-no-relaxed";
  string fixedvars = "-fixed";
  string onlybreakers = "-print-only-breakers";
  string generatorfile = "-store-sym";
  string symmetryinput = "-addsym";
  string aspinput = "-asp";
  string pbinput = "-pb";
  string xnfinput = "-xnf";
  string logfile = "-logfile";
}

void printUsage() {
  std::clog << "BreakID version " << VERSION << std::endl;
  std::clog << "Usage: ./BreakID <cnf-file> " <<
          "[" << options::help << "] " <<
          "[" << options::instancefile << " <file with instance>]" <<
          "[" << options::nointch << "] " <<
          "[" << options::nobinary << "] " <<
          "[" << options::nosmall << "] " <<
          "[" << options::norelaxed << "] " <<
          "[" << options::formlength << " <nat>] " <<
          "[" << options::timelim << " <nat>] " <<
          "[" << options::verbosity << " <nat>] " <<
          "[" << options::fixedvars << " <file with fixed vars>] " <<
          "[" << options::onlybreakers << "] " <<
          "[" << options::generatorfile << " <path to store symmetry generators>] " <<
          "[" << options::symmetryinput << "<file with symmetry info>] " <<
          "[" << options::aspinput << "] " <<
          "[" << options::pbinput << " <nat>] " <<
          "[" << options::xnfinput << "] " <<
          "[" << options::logfile << " <file to store veriPB log in>]  " <<
          "\n";
  std::clog << "\nOptions:\n";
  std::clog << options::help << "\n  ";
  std::clog << "Display this help message instead of running BreakID.\n";
  std::clog << options::instancefile << "\n";
  std::clog << "Read instance from a file instead of input stream.\n";
  std::clog << options::nointch << "\n  ";
  std::clog << "Disable detection and breaking of row interchangeability.\n";
  std::clog << options::nobinary << "\n  ";
  std::clog << "Disable construction of additional binary symmetry breaking clauses based on stabilizer subgroups.\n";
  std::clog << options::nosmall << "\n  ";
  std::clog << "Disable compact symmetry breaking encoding, use Shatter's encoding instead.\n";
  std::clog << options::norelaxed << "\n  ";
  std::clog << "Disable relaxing constraints on auxiliary encoding variables, use longer encoding instead.\n";
  std::clog << options::formlength << " <default: " << symBreakingFormLength << ">\n  ";
  std::clog << "Limit the size of the constructed symmetry breaking formula's, measured as the number of auxiliary variables introduced. <-1> means no symmetry breaking.\n";
  std::clog << options::timelim << " <default: " << timeLim << ">\n  ";
  std::clog << "Upper limit on time spent by Saucy detecting symmetry measured in seconds.\n";
  std::clog << options::verbosity << " <default: " << verbosity << ">\n  ";
  std::clog << "Verbosity of the output. <0> means no output other than the CNF augmented with symmetry breaking clauses.\n";
  std::clog << options::fixedvars << " <default: none>\n  ";
  std::clog << "File with a list of variables that should be fixed, separated by whitespace.\n";
  std::clog << options::onlybreakers << "\n  ";
  std::clog << "Do not print original theory, only the symmetry breaking clauses.\n";
  std::clog << options::generatorfile << "\n  ";
  std::clog << "Store the detected symmetry generators in the given file.\n";
  std::clog << options::symmetryinput << " <default: none>\n  ";
  std::clog << "Pass a file with symmetry generators or row-interchangeable matrices to use as additional symmetry information. Same format as BreakID's output by "
            << options::generatorfile << ".\n";
  std::clog << options::aspinput << "\n  ";
  std::clog << "Parse input in the smodels-lparse intermediate format instead of DIMACS.\n";
  std::clog << options::pbinput << " <default:unset>\n  ";
  std::clog
          << "When this variable is set, we assume input in the OPB format. The value of the integer indicates the number of variables grouped together in one PB lex-leader constraint. In particular: 1 means almost clausal encoding and infinity means one PB constraint with exponentially sized coefficients. Special case is 0, which takes the PB encoding of symmetry breaking CNF clauses.\n";
  std::clog << options::xnfinput << "\n  ";
  std::clog << "When this variable is set, we assume input in the XNF format.\n";
  std::clog << options::logfile << "\n  ";
  std::clog << "Store a log in veriPB format in the provided file <default:unset>.\n";
  gracefulError("");
}

void parseOptions(int argc, char *argv[]) {
  for (int i = 1; i < argc; ++i) {
    string input = argv[i];
    if (0 == input.compare(options::nobinary)) {
      useBinaryClauses = false;
    } else if (0 == input.compare(options::nointch)) {
      useMatrixDetection = false;
    } else if (0 == input.compare(options::onlybreakers)) {
      onlyPrintBreakers = true;
    } else if (0 == input.compare(options::generatorfile)) {
      ++i;
      generatorFile = argv[i];
    } else if (0 == input.compare(options::logfile)) {
      ++i;
      logging = true;
      logfile = argv[i];
      logstream.open(logfile, ios::out);
      if (logstream.fail()) {
        gracefulError("Could not write to "+logfile);
      }
    } else if (0 == input.compare(options::nosmall)) {
      useShatterTranslation = true;
    } else if (0 == input.compare(options::norelaxed)) {
      useFullTranslation = true;
    } else if (0 == input.compare(options::formlength)) {
      ++i;
      symBreakingFormLength = stoi(argv[i]);
    } else if (0 == input.compare(options::timelim)) {
      ++i;
      timeLim = stoi(argv[i]);
    } else if (0 == input.compare(options::verbosity)) {
      ++i;
      verbosity = stoi(argv[i]);
    } else if (0 == input.compare(options::help)) {
      printUsage();
    } else if (0 == input.compare(options::fixedvars)) {
      ++i;
      string filename = argv[i];
      setFixedLits(filename);
    } else if (0 == input.compare(options::aspinput)) {
      aspinput = true;
    } else if (0 == input.compare(options::pbinput)) {
      ++i;
      nbPBConstraintsToGroup = stoi(argv[i]);
    } else if (0 == input.compare(options::xnfinput)) {
      xnfinput = true;
    } else if (0 == input.compare(options::symmetryinput)) {
      ++i;
      inputSymFile = argv[i];
    } else if (0 == input.compare(options::instancefile)) {
      ++i;
      inputFile = argv[i];
    }
  }

  if (verbosity > 1) {
    std::clog << "Options used: " <<
            options::formlength << " " << symBreakingFormLength << " " <<
            options::timelim << " " << timeLim << " " <<
            options::verbosity << " " << verbosity << " " <<
            (useMatrixDetection ? "" : options::nointch) << " " <<
            (useBinaryClauses ? "" : options::nobinary) << " " <<
            (onlyPrintBreakers ? options::onlybreakers : "") << " " <<
            (generatorFile != "" ? options::generatorfile  + " " + generatorFile : " ") << " " <<
            (logfile != "" ? options::logfile + " " + logfile : " ") << " " <<
            (inputSymFile != "" ? options::symmetryinput + " " + inputSymFile : " ") << " " <<
            (inputFile != "" ? options::instancefile + " " + inputFile : " ") << " " <<
            (useShatterTranslation ? options::nosmall : "") << " " <<
            (useFullTranslation ? options::norelaxed : "") << " " <<
            (fixedLits.size()==0 ? "" : options::fixedvars) << " " << //TODO PRINT AND STORE FIXED LITS FILE NAME
            (aspinput ? options::aspinput : " ") << " " <<
            (xnfinput ? options::xnfinput : " ") << " " <<
            std::endl;
  }
}

// ==== main

void checkOptionsCompatible(){
  if(logging && useBinaryClauses){
    std::stringstream ss;
    ss << "Please disable binary clause detection when using proof logging. This is not implemented yet. It can be done with " << options::nobinary;
    gracefulError(ss.str());
  }
  if((aspinput | (nbPBConstraintsToGroup != -1 ) ) && logging){
     std::stringstream ss;
    ss << "Logging is only implemented with CNF input (not with asp or pb input yet)";
    gracefulError(ss.str());
  }

  if(aspinput && onlyPrintBreakers){
    std::stringstream ss;
    ss << "Options " << options::aspinput << " and " << options::onlybreakers << " are incompatible since asp output is intertwined with extra clauses.\n";
    gracefulError(ss.str());
  }
  if (xnfinput && (nbPBConstraintsToGroup != -1 || aspinput)) {
    std::stringstream ss;
    ss << "Option " << options::xnfinput << " is incompatible with options " << options::aspinput << " and " << options::pbinput << ".\n";
    gracefulError(ss.str());
  }
}

sptr<Specification> readTheory(std::istream& input){
    sptr<Specification> theory;
    if(aspinput){
        theory = make_shared<LogicProgram>(input);
    } else if (nbPBConstraintsToGroup != -1) {
        theory = make_shared<PB>(input);
    } else if (xnfinput) {
        theory = make_shared<XNF>(input);
    } else {
        theory = make_shared<CNF>(input);
    }
    return theory;
}

int main(int argc, char *argv[]) {
  time(&startTime);

  parseOptions(argc, argv);
  checkOptionsCompatible();

  sptr<Specification> theory;
  if(inputFile!=""){
      if (verbosity > 0) {
          std::clog << "*** Reading instance: " << inputFile << std::endl;
      }
      std::ifstream input(inputFile);
      if (input.fail()) {
          gracefulError("No instance file found.");
      }
      theory = readTheory(input);
  }else{
      if (verbosity > 0) {
          std::clog << "*** Reading instance from input stream." << std::endl;
      }
      theory = readTheory(std::cin);
  }

  if (verbosity > 3) {
    theory->getGraph()->print();
  }

  if (verbosity > 0) {
    std::clog << "**** symmetry generators detected: " << theory->getGroup()->getSize() << std::endl;
    if (verbosity > 2) {
      theory->getGroup()->print(std::clog);
    }
  }
  
  addInputSym(theory->getGroup());
  
  if(verbosity>0){
    std::clog << "*** Detecting subgroups..." << std::endl;
  }
  vector<sptr<Group> > subgroups;
  theory->getGroup()->getDisjointGenerators(subgroups);
  if (verbosity > 0) {
    std::clog << "**** subgroups detected: " << subgroups.size() << std::endl;
  }

  if (verbosity > 1) {
    for (auto grp : subgroups) {
      std::clog << "group size: " << grp->getSize() << " support: " << grp->getSupportSize() << std::endl;
      if (verbosity > 2) {
        grp->print(std::clog);
      }
    }
  }
  
  theory->cleanUp(); // improve some memory overhead
  
  unsigned int totalNbMatrices = 0;
  unsigned int totalNbRowSwaps = 0;
  if(logging){
    nbConstraintsInProofLog = theory->getOriginalSize(); 
    logstream << "pseudo-Boolean proof version 1.2\n"
     << "* load formula\n"
     << "f " << nbConstraintsInProofLog << "\n\n";

  }
  
  Breaker brkr(theory);
  for (auto grp : subgroups) {
    if (grp->getSize() > 1 && useMatrixDetection) {
      if (verbosity > 1) {
        std::clog << "*** Detecting row interchangeability..." << std::endl;
      }
      theory->setSubTheory(grp);
      grp->addMatrices();
      totalNbMatrices += grp->getNbMatrices();
      totalNbRowSwaps += grp->getNbRowSwaps();
    }
    if(symBreakingFormLength > -1){
        if (verbosity > 1) {
            std::clog << "*** Constructing symmetry breaking formula..." << std::endl;
        }
        grp->addBreakingClausesTo(brkr);  
    }// else no symmetry breaking formulas are needed :)
    grp.reset();
  }
  
  if (verbosity > 0) {
    std::clog << "**** matrices detected: " << totalNbMatrices << std::endl;
    std::clog << "**** row swaps detected: " << totalNbRowSwaps << std::endl;
    std::clog << "**** extra binary symmetry breaking clauses added: " << brkr.getNbBinClauses() << "\n";
    std::clog << "**** regular symmetry breaking clauses added: " << brkr.getNbRegClauses() << "\n";
    std::clog << "**** row interchangeability breaking clauses added: " << brkr.getNbRowClauses() << "\n";
    std::clog << "**** total symmetry breaking clauses added: " << brkr.getAddedNbClauses() << "\n";
    std::clog << "**** auxiliary variables introduced: " << brkr.getAuxiliaryNbVars() << "\n";
    std::clog << "*** Printing resulting theory with symmetry breaking clauses..." << std::endl;
  }

  brkr.print();

  if(generatorFile!=""){
    if (verbosity > 0) {
      std::clog << "*** Printing generators to file "+generatorFile << std::endl;
    }
    ofstream fp_out;
    fp_out.open(generatorFile, ios::out);
    if (fp_out.fail()) {
      gracefulError("Could not write to "+generatorFile);
    }

    for (auto grp : subgroups) {
      grp->print(fp_out);
    }
    fp_out.close();
  }
  if(logging){
    logstream.close();
  }

  return 0;
}
