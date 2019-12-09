// #include "TH3.h"
// #include "TH2.h"
// #include "TH1.h"
// #include "TFile.h"
// #include <TROOT.h>
// #include "TTree.h"
// #include "TProfile.h"
// #include "TCanvas.h"
// #include "TFitResultPtr.h"
// #include "TMath.h"
#include "dataTypes.h"
#include<algorithm>
#include<iostream>
#include<iomanip>
#include<math.h>
#include<fstream>
#include<sstream>
using namespace std;

int binomialCoefficients(int l, int k) {
   if (k == 0 || k == l)
   return 1;
   return binomialCoefficients(l - 1, k - 1) + binomialCoefficients(l - 1, k);
}
double legPoly(int l, double x)
{
  double pl = 0.;
  for (int k=0;k<l+1;k++)
  {
    pl+=pow(binomialCoefficients(l,k),2)*pow(x-1,l-k)*pow(x+1,k);
  }
  return pl/pow(2,l);
}
uint16_t riplLevel::getDecayLevel(double a_randNum)
{
  // for (double ee : m_branchRatio) cout << "Br = " << ee << ", ";
  // cout << m_branchRatio.size() <<"\n";
  // cout << "\n";
  int nextLevelIndex = distance(m_branchRatio.begin(),upper_bound(m_branchRatio.begin(),m_branchRatio.end(),a_randNum));
  // cout << "next level = " << nextLevelIndex<<", " << m_nextLevelID[nextLevelIndex]<<"\n";
  return m_nextLevelID[nextLevelIndex];
}
void riplLevel::prepareLevel()
{
  vector<double> branchSort;
  double BR_tot = 0;
  for (double BR : m_branchRatio)
  {
    BR_tot+=BR;
  }
  double BR_cumul=0.;
  for (double BR : m_branchRatio)
  {
    branchSort.push_back(BR/BR_tot + BR_cumul);
    BR_cumul+=BR/BR_tot;
  }
  m_branchRatio = branchSort;
}
void RIPL::readRIPL(double BR_cutoff)
{
  ifstream fRIPL ("../fe56_ripl.txt");
  string line;
  vector<vector<double>> outp;
  int ensdfLevel = 1;
  if (fRIPL.is_open())
  {
    for (int i=0; i<13; i++) {
      getline (fRIPL, line);
      // if (i==2) setNumLevels(stoi(line.substr(20,50)));
      // if (i==2) numLvls = atoi();
      // else if (i==3) numGs = atoi();
    }
    double EL = -1.;
    uint16_t EL_ID = -1;
    double EG = -1.;
    double PE = -1.;
    double PG = -1.;
    int numGammas = -1;
    vector<vector<double>> EG_EL;
    // vector<double> BR;
    bool endOfData=false;
    vector<double> PG_lvl;
    double PG_sum = 0.;
    vector<uint16_t> levelIDs;
    while (!endOfData) 
    {
      /// new level
      getline (fRIPL, line);
      // cout << "LINE " << line << "\n";
      if (line[0]=='-') endOfData=true;
      else if (line[2]!=' ') 
      {
        // cout << EL << ": " << PG_sum << "\n";
        if (numGammas>0 && PG_sum>0)
        {
          m_data[EL_ID].setEnergy(EL);
          m_data[EL_ID].setENSDFLevel(ensdfLevel);
          m_data[EL_ID].setLevelID(levelIDs);
          m_data[EL_ID].setBR(PG_lvl);
          PG_lvl.clear();
          levelIDs.clear();
          ensdfLevel++;
          PG_sum=0;
        }
        EL = stof(line.substr(3,11));
        
        EL_ID = (uint16_t)stoi(line.substr(0,3));
        numGammas = stoi(line.substr(36,2));
        // cout << EL_ID << ": " << numGammas<<"\n";


      }
      else if (line[0]!='-') 
      {
        EG = stof(line.substr(43,11));
        PG = stof(line.substr(54,11));

        PE = stof(line.substr(65,11));

        // m_data[EL_ID].m_nextLevelID.push_back(atoi(line.substr(0,43).c_str()));
        // m_data[EL_ID].m_branchRatio.push_back(PG);
        if (numGammas>0 )
        {
          // m_data[EL_ID].addLevelID((uint16_t)stoi(line.substr(0,43)));
          // m_data[EL_ID].addBR(PG);
          PG_lvl.push_back(PG);
          PG_sum+=PG;
          levelIDs.push_back((uint16_t)stoi(line.substr(0,43)));
        }
      }
    }
  }
  else cout << "could not open RIPL data\n";
  cout << "Number of RIPL Levels processed = " <<ensdfLevel << "\n";
  setNumLevels(ensdfLevel);
  for (auto& pair : m_data)
  {
    pair.second.prepareLevel();
    // cout << pair.first << ", " << pair.second.getENSDFLevel() << ", " << pair.second.getEnergy()<<"\n";
  }
}

vector<double> RIPL::doCascade(uint16_t a_initLevel, randNum* a_randGen)
{
  vector<double> gammas;
  uint16_t nextLevel =0;
  uint16_t prevLevel = a_initLevel;
  while (nextLevel!=1)
  // for (int i=0;i<5;i++)
  {
    nextLevel = m_data[prevLevel].getDecayLevel(a_randGen->next());
    // cout << prevLevel << "--> " << nextLevel << "\n";
    gammas.push_back(m_data[prevLevel].getEnergy()-m_data[nextLevel].getEnergy());
    prevLevel = nextLevel;
  }
  return gammas;
}
uint16_t RIPL::getLevelNumber(string a_level, randNum* a_randGen)
{
  uint16_t level = stoi(a_level)-49;
  int level_ensdf = stoi(a_level)-50;

  for (auto& pair : m_data)
  {
    if (level_ensdf==pair.second.getENSDFLevel()) level = pair.first;
  }

  
  if (level>40)
  {
    level = (uint16_t)floor((m_numLevels-40)*a_randGen->next() + 40);
    // cout << "begin: "<< level << ", ";
    std::map<uint16_t, riplLevel>::iterator it = m_data.find(level);
    while (it==m_data.end())
    {
      // cout << level << ", " ;
      level = (uint16_t)floor((m_numLevels-40)*a_randGen->next() + 40);
      it = m_data.find(level);
    }
  }
  // else if (it==m_data.end())
  // cout << "\n";
  return level;
}
double RIPL::getLevelEnergy(uint16_t a_level)
{
  return m_data[a_level].getEnergy();
}

double angDistLeg::sample(double a_rand)
{
  vector<double> mu;
  vector<double> prob;
  double prob_tot=0.;
  for (int p=0; p<200; p++)
  {
    mu.push_back(-1.+2.*p/200);
    double pl =0.;
    for (int l =0; l<m_coeffs.size(); l++)
    {
      pl += m_coeffs[l]*legPoly(l,-1.+2.*p/200);
    }
    pl+=1.;
    prob_tot+=pl;
    prob.push_back(pl);
  }
  double norm_prob_sum=0.;
  for (int iP=0; iP<200; iP++)
  {
    double norm_prob = prob[iP]/prob_tot;
    double prob_cur  = norm_prob+norm_prob_sum;
    norm_prob_sum   += norm_prob;
    if (a_rand<prob_cur) return mu[iP];
  }
}

void angDistLeg::printData()
{
  for (double coeff: m_coeffs) cout << coeff << ", ";
  cout << "\n";
}