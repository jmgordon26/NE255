// project includes
#include "readXSec.h"
readXSec::readXSec()
{
  ;
}

void readXSec::initPartialInel(string a_isotope)
{
  string fname = "../FeXsec/"+a_isotope+"_InelF01";
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
  int mf =3;
  // cout << "first MT = " << mt << "\n";
  /// loop over level cross sections (MF=3)
  while(true)
  {
    if (mf==3)
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
      mf = stoi(s);
      // if (stoi(s)!=3) break;
      getline(ifs, s);
      lineNum++;
      mt = stoi(s.substr(2,2));    
    }
    else if (mf==4)
    {
      while (s.size()<20) getline(ifs, s);
      // vector<double> coeffs;
      double energy = 0.;
      angDistLeg* coeffs = new angDistLeg();
      while (s.size()>3)
      {
        if (stof(s.substr(0,12))==0)
        {
          if (coeffs->m_coeffs.size()>0) 
          {
            // cout << "\n"<<to_string(mt)<<": NEW ENERGY "<<energy<<"\n";
            // coeffs->printData();
            m_partialInelAng[to_string(mt)][energy] = (*coeffs);
            coeffs->m_coeffs.clear();
          }
          energy = (double)stof(s.substr(14,12));
        }
        else
        {
          for (int i=0; i<s.size()/14; i++) 
          {
            coeffs->m_coeffs.push_back((double)stof(s.substr(i*14,14)));
            // cout << (double)stof(s.substr(i*14,14)) << ", ";
          }
        }
        getline(ifs,s);
        lineNum++;
      }
      // m_partialInelAng[to_string(mt)]
      getline(ifs, s);
      lineNum++;
      // check for MF number
      mf = stoi(s);
      // if (stoi(s)!=3) break;
      getline(ifs, s);
      lineNum++;
      mt = stoi(s.substr(2,2));    
    }
    else break;
  }
  ifs.close();
  // for (auto& pair : m_partialInelAng["89"])
  // {
  //   cout << pair.first << ": " << pair.second.m_coeffs.size() << "\n";
  // }
}

void readXSec::initTotal(string a_isotope)
{
  string fname;
  map<string, string> fnames {{"el","_ElXsec"}, {"inel","_InelXsec"}, {"cap","_CapXsec"}};//, {"fiss","_FisXsec"} };

  string s;
  for (auto& pair : fnames)
  {
    fname = "../FeXsec/"+a_isotope+pair.second;
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
      // cout << pair.first << ": " << m_total[pair.first].energy.size() <<"\n";
    }
    ifs.close();
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
  string fname = "../FeXsec/"+a_isotope+"_CapProd";
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
  ifs.close();
}

double readXSec::pickInelAngle(string a_level, double a_incidentEnergy, double a_rand)
{
  // cout << a_level << "\n";
  auto iter = m_partialInelAng.find(a_level);
  if (iter!=m_partialInelAng.end())
  {
    auto iter = m_partialInelAng[a_level].begin();
    double a_dist_cur = abs((*iter).first - a_incidentEnergy);
    
    double a_dist_prev = abs((*iter).first - a_incidentEnergy);
    // cout << "BeamEnergy = " << a_incidentEnergy << ", " <<(*iter).first<< "\n";
    while (a_dist_prev>=a_dist_cur)
    {
      a_dist_prev = a_dist_cur;
      iter = next(iter,1);
      a_dist_cur = abs((*iter).first - a_incidentEnergy);
      // cout << (*iter).first << ", " << abs((*iter).first - a_incidentEnergy) << "\n";
    }
    iter = prev(iter,1);
    // cout << "found energy = " << (*iter).first << "\n";
    double mu = (*iter).second.sample(a_rand);
    return mu;    
  }
  else 
  {
    return a_rand*2.-1;
  }
}

void readXSec::readElasticAngle(string a_isotope)
{
  string fname = "../FeXsec/"+a_isotope+"_ElFS";
  ifstream ifs {fname.c_str()};
  string s;
  while (s.size()<20) getline(ifs, s);
  // vector<double> coeffs;
  double energy = 0.;
  angDistLeg* coeffs = new angDistLeg();
  while (s.size()>20)
  {
    if (stof(s.substr(0,12))==0)
    {
      if (coeffs->m_coeffs.size()>0) 
      {
        // cout << "\nNEW ENERGY "<<energy<<"\n";
        // coeffs->printData();
        m_elasticAng[energy] = (*coeffs);
        m_elasticAngEnergy.push_back(energy);
        coeffs->m_coeffs.clear();
      }
      energy = (double)stof(s.substr(14,12));
    }
    else
    {
      for (int i=0; i<s.size()/14; i++) 
      {
        coeffs->m_coeffs.push_back((double)stof(s.substr(i*14,14)));
        // cout << (double)stof(s.substr(i*14,14)) << ", ";
      }
    }
    getline(ifs,s);

  }
  cout <<"elastic angle size = "<< m_elasticAng.size() << "\n";
  ifs.close();
}

double readXSec::pickElasticAngle(double a_incidentEnergy, double a_rand)
{
  // cout << a_level << "\n";
  if (a_incidentEnergy>1e6)
  {
    // auto iter = m_elasticAng.end();
    // iter = prev(iter,1);
    // double a_dist_cur = abs((*iter).first - a_incidentEnergy);

    // double a_dist_prev = abs((*iter).first - a_incidentEnergy);
    // // cout << "BeamEnergy = " << a_incidentEnergy << ", " <<(*iter).first<< "\n";
    // while (a_dist_prev>=a_dist_cur)
    // {
    //   a_dist_prev = a_dist_cur;
    //   iter = prev(iter,1);
    //   a_dist_cur = abs((*iter).first - a_incidentEnergy);
    //   // cout << (*iter).first << ", " << abs((*iter).first - a_incidentEnergy) << "\n";
    // }
    // iter = next(iter,1);
    // cout << "found energy = " << (*iter).first << "\n";
    double energy = (*upper_bound(m_elasticAngEnergy.begin(),m_elasticAngEnergy.end(),a_incidentEnergy));
    // double mu = (*iter).second.sample(a_rand);
    double mu = m_elasticAng[energy].sample(a_rand);
    return mu;    
  }
  else
  {
    // auto iter = m_elasticAng.begin();
    // iter = next(iter,1);
    // double a_dist_cur = abs((*iter).first - a_incidentEnergy);

    // double a_dist_prev = abs((*iter).first - a_incidentEnergy);
    // // cout << "BeamEnergy = " << a_incidentEnergy << ", " <<(*iter).first<< "\n";
    // while (a_dist_prev>=a_dist_cur)
    // {
    //   a_dist_prev = a_dist_cur;
    //   iter = next(iter,1);
    //   a_dist_cur = abs((*iter).first - a_incidentEnergy);
    //   // cout << (*iter).first << ", " << abs((*iter).first - a_incidentEnergy) << "\n";
    // }
    // iter = prev(iter,1);
    // cout << "found energy = " << (*iter).first << "\n";
    double energy = (*upper_bound(m_elasticAngEnergy.begin(),m_elasticAngEnergy.end(),a_incidentEnergy));
    // double mu = (*iter).second.sample(a_rand);
    double mu = m_elasticAng[energy].sample(a_rand);
    return mu;     
  }

}
vector<double> readXSec::readENSDF(string a_isotope, double a_rand)
{

  string fname = "../../FeXsec/ensdf_"+a_isotope.substr(0,a_isotope.find(" "))+".txt";
  // cout << fname << "\n";
  ifstream ifs {fname.c_str()};
  if (!ifs) 
  {
    // cout << fname<<"\n";
    return {};
  }
  string line;
  getline(ifs, line);
  string daughterNuc = line.substr(1,6)+"G";
  string daughterNucLvl = line.substr(1,6)+"L";
  double tot = 0.;
  vector<double> gammas_out;
  vector<double> gammas;
  vector<double> ratios;
  
  while (line.size()>1)
  {
    if (line.substr(1,7)==daughterNucLvl)
    {
      double level = stof(line.substr(9,10));
      // cout <<"level = "<< level << "\n";
    }
    else if (line.substr(1,7)==daughterNuc)
    {
      double energy = stof(line.substr(9,10));
      // cout << energy << "\n";
      double br = stof(line.substr(21,5));

      // cout << br << "\n";
      // tot+=br;
      if (br>=100) gammas_out.push_back(energy);
      else 
      {
        tot+=br;
        gammas.push_back(energy);
        ratios.push_back(br);
      }
    }
    getline(ifs,line);
  }
  // cout <<"total = "<< tot << "\n";
  double tot_cum = 0.;
  for (int iG =0; iG<gammas.size(); iG++)
  {
    double br = ratios[iG]/tot;

    tot_cum+=br;
    // cout << br << ", " <<tot_cum << ", " <<a_rand<<"\n";
    if (tot_cum>a_rand) {gammas_out.push_back(gammas[iG]); break;}
  }
  // 
  return gammas_out;
}