#include "randNum.h"
#include <vector>
#include <map>
#include <string>

/// class to hold cross section data
/// in energy vs. xsec
class xSec
{
public:
  xSec() {;}
  std::vector<double> energy;
  std::vector<double> xsec;

};
/// classes to hold angular distribution data
/// will have histogrammed and legendre data
class angDist
{
public:
  angDist() {;}
  std::vector<double> m_coeffs;
protected:
  double sample(double a_rand);

};
class angDistLeg : public angDist
{
public:
  angDistLeg() {;}
  std::vector<double> m_coeffs;

  double sample(double a_rand);
  void printData();
};
class angDistHist : public angDist
{
public:
  angDistHist() {;}
  std::vector<double> m_mu;
  std::vector<double> m_prob;

  //double sample(double a_rand);
};
/// represents decay gammas from a given level
class riplLevel
{
public:
  /// default constructor
  riplLevel() {;}
  /// fix branching ratios for random sampling
  void prepareLevel();
  /// returns level ID of next level - randomly samp'd BR's
  uint16_t getDecayLevel(double a_randNum);
  /// set the energy of the level [MeV]
  void setEnergy(double a_energy) {m_energy = a_energy;}
  /// set the RIPL level number
  void setENSDFLevel(int a_ensdfLevel) {m_ensdfLevel = a_ensdfLevel;}
  int getENSDFLevel() {return m_ensdfLevel;}
  /// returns energy of the level [MeV]
  double getEnergy() {return m_energy;};
  /// add level number of decay path
  // void addLevelID(uint16_t a_levelIDs) {m_nextLevelID.push_back(a_levelIDs);}
  void setLevelID(vector<uint16_t> a_levelIDS) {m_nextLevelID = a_levelIDS;}
  /// add branching ratio of decay gamma
  // void addBR(double a_BR) {m_branchRatio.push_back(a_BR);}
  void setBR(vector<double> a_BR) {m_branchRatio = a_BR;}

private:
  int m_ensdfLevel = 0;
  double m_energy=0;
  std::vector<uint16_t> m_nextLevelID;
  std::vector<double> m_branchRatio;
};

/// holds RIPL data for an isotope
/// also does cascade
class RIPL
{
public:
  /// default constructor
  RIPL() {;}
  /// read data from file
  void readRIPL(double BR_cutoff);
  /// returns energy of level a_level
  double getLevelEnergy(uint16_t a_level);
  /// return RIPL level number
  uint16_t getLevelNumber(std::string a_level, randNum* a_randGen);
  /// set number of levels read
  void setNumLevels(int a_numLevels) {m_numLevels = a_numLevels;}
  /// does cascade - returns gammas produced [MeV]
  std::vector<double> doCascade(uint16_t a_initLevel, randNum* a_randGen);
private:
  /// level data
  std::map<uint16_t, riplLevel> m_data;
  /// number of levels
  int m_numLevels =0;
};