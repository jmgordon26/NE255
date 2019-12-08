// project includes
// #include "randNum.h"
#include "dataTypes.h"
// c++ includes
#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<cmath>
#include<string>
#include<vector>
#include<map>
#include<algorithm>
// root includes
// #include "TGraph.h"

using namespace std;

// simple G4ENDL file reader

class readXSec {

public:
  readXSec();

  void setRandGenerator(randNum* a_randGen) {m_rand = a_randGen;}

  void initPartialInel(string a_isotope);
  void initTotal(string a_isotope);
  void readCapture(string a_isotope);

  string pickReaction(double a_energy);

  string pickInelLevel(double a_energy);

private:
  randNum* m_rand;
  /// map of total cross sections
  map<string, xSec> m_total;
  /// map of lineNumber of partial inelastic
  map<int, vector<uint16_t>> m_partialInel;
  map<string, xSec> m_partialInelx;
  map<string, xSec> m_isotopeProd;

};