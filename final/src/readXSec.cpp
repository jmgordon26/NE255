// project includes
#include "readXSec.h"
readXSec::readXSec()
{
  ;
}

void readXSec::initPartialInel(string a_isotope)
{
  string fname = "../../FeXsec/"+a_isotope+"_InelF01";
  ifstream ifs {fname.c_str()};
  map<int, vector<vector<double>>> mt_map;
  // ifstream ifs {"~/ENDF-VIII.0/Inelastic/F01/26_56_Iron"};
  if (!ifs) {cout << "couldn't open file :(\n"; return ;}
  // else cout << "file opened successfully!\n";
  string s;
  uint16_t lineNum=0;
  getline(ifs, s);
  lineNum++;
  while (s.size()==3) 
  {
    getline(ifs, s);
    lineNum++;
  }

  int mt = stoi(s.substr(2,2));
  // cout << "first MT = " << mt << "\n";
  /// loop over level cross sections (MF=3)
  while(true)
  {
    getline(ifs, s);
    lineNum++;
    getline(ifs, s);
    lineNum++;
    vector<double> energy;
    vector<double> xsec;
    uint16_t lineNumStart = lineNum;
    while (s.size()>3)
    {
      for (int iD=0;iD<s.size()/28;iD++)
      {
        energy.push_back((double)stof(s.substr(2+28*iD,12)));
        xsec.push_back((double)stof(s.substr(16+28*iD,12)));        
      }

      getline(ifs, s);
      lineNum++;
    }
    m_partialInel[mt] = {lineNumStart,lineNum};
    // mt_map[mt] = {energy,xsec};
    if (mt>4)
    {
      m_partialInelx[to_string(mt)].energy = energy;
      m_partialInelx[to_string(mt)].xsec = xsec;    
    }
      getline(ifs, s);
    lineNum++;
    // check for MF number
    if (stoi(s)!=3) break;
    getline(ifs, s);
    lineNum++;
    mt = stoi(s.substr(2,2));
  }

}

void readXSec::initTotal(string a_isotope)
{
  string fname;
  map<string, string> fnames {{"el","_ElXsec"}, {"inel","_InelXsec"}, {"cap","_CapXsec"}};//, {"fiss","_FisXsec"} };

  string s;
  for (auto& pair : fnames)
  {
    fname = "../../FeXsec/"+a_isotope+pair.second;
    cout << "working on " << pair.second <<"\n";
    ifstream ifs {fname.c_str()};
    if (!ifs) {cout << "couldn't open files :("; return;}
    getline(ifs, s);
    while (s.size()<10) getline(ifs, s);
    /// loop over entries in file
    while (s.size()>10)
    {
      vector<double> energy;
      vector<double> xsec;
      while (s.size()>3)
      {
        for (int iD=0;iD<s.size()/28;iD++)
        {
          energy.push_back((double)stof(s.substr(2+28*iD,12)));
          xsec.push_back((double)stof(s.substr(16+28*iD,12)));        
        }
        getline(ifs, s);
      }
      m_total[pair.first].energy = energy;
      m_total[pair.first].xsec = xsec;
      cout << pair.first << ": " << m_total[pair.first].energy.size() <<"\n";
    }
  }
}

string readXSec::pickReaction(double a_energy)
{
  double rand1 = m_rand->next();
  map<string, double> crossSections;
  double totalXsec=0.;
  // cout << rand1 << "\n";
  for (auto& pair : m_total)
  {
    vector<double> energy = pair.second.energy;
    int index = distance(energy.begin(),upper_bound(energy.begin(),energy.end(), a_energy));
    if (index>0)
    {
      double xsec = pair.second.xsec[index-1];
      crossSections[pair.first]=xsec;
      totalXsec+=xsec;
      // cout << pair.first << ": " << xsec << ": " << energy[index] << ", ";  
    }

    // if (pair.first=="cap") cout << "total = " << xsec;
    
  }
  double totalXsec_cap=0.;
  for (auto& pair : m_isotopeProd)
  {
    vector<double> energy = pair.second.energy;
    int index = distance(energy.begin(),upper_bound(energy.begin(),energy.end(), a_energy));
    if (index>0)
    {
      double xsec = pair.second.xsec[index-1];
      crossSections[pair.first]=xsec;
      totalXsec+=xsec;
      // cout << pair.first << ": " << xsec << ": " << energy[index] << ", ";  
    }
  }
  // crossSections["iso"]=totalXsec_cap;
  // totalXsec+=totalXsec_cap;
  // cout << "\n";
  // cout << " and capIsotope = " << totalXsec_cap << "\n";
  double normXsecSum=0;
  for (auto& pair : crossSections)
  {
    double normXsec = pair.second/totalXsec;
    pair.second=normXsec+normXsecSum;
    normXsecSum+=normXsec;
    // cout << "\t"<<pair.first << ": "<< pair.second << "\n";
    if (rand1<pair.second) return pair.first;
  }
}

string readXSec::pickInelLevel(double a_energy)
{
  double rand1 = m_rand->next();
  map<string, double> crossSections;
  double totalXsec=0.;
  // cout << rand1 << "\n";
  for (auto& pair : m_partialInelx)
  {
    vector<double> energy = pair.second.energy;
    int index = distance(energy.begin(),upper_bound(energy.begin(),energy.end(), a_energy));
    double xsec = pair.second.xsec[index];
    crossSections[pair.first]=xsec;
    totalXsec+=xsec;
    // cout << pair.first << ": " << xsec << ", ";
  }
  // cout << "\n";
  double normXsecSum=0;
  for (auto& pair : crossSections)
  {
    double normXsec = pair.second/totalXsec;
    pair.second=normXsec+normXsecSum;
    normXsecSum+=normXsec;
    // cout << "\t"<<pair.first << ": "<< pair.second << "\n";
    if (rand1<pair.second) return pair.first;
  }
}

void readXSec::readCapture(string a_isotope)
{
  string fname = "../../FeXsec/"+a_isotope+"_CapProd";
  ifstream ifs {fname.c_str()};

  string isotope="";
  vector<double> energy;
  vector<double> xsec;

  if (!ifs) {cout << "couldn't open file :(\n"; return ;}

  string line;
  while (getline(ifs,line))
  {
    if (line[0]!=' ') 
    {
      if (xsec.size()>0)
      {
        m_isotopeProd[isotope].xsec = xsec;
        m_isotopeProd[isotope].energy = energy;
        energy.clear();
        xsec.clear();        
      }

      isotope=line;
    }
    else if (line.size()>20)
    {
      for (int iD=0;iD<line.size()/22;iD++)
      {
        energy.push_back((double)stof(line.substr(2*iD*11,11)));
        xsec.push_back((double)stof(line.substr((2*iD+1)*11,11)));
      }
    }
  }
}

// void pickC