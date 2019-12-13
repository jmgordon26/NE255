// #include "dataTypes.h"
// #include "randNum.h"
#include "readXSec.h"
// #include "extraMath.h"

#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<cmath>
#include<string>
#include <cstring>
#include<vector>
#include<map>
#include<algorithm>
#include <stdlib.h>
using namespace std;

std::vector<double> initPosition(double a_rad, randNum* a_randGen)
{
  double r1 = a_rad*(2*a_randGen->next()-1);
  double r2 = a_rad*(2*a_randGen->next()-1);

  double radius = sqrt(r1*r1+r2*r2);
  while (radius > a_rad)
  {
    r2 = a_rad*(2*a_randGen->next()-1);
    radius = sqrt(r1*r1+r2*r2);
  }
  return {r1,r2};
}

int main(int argc, char *argv[])
{
  /// process args: numParticles, verbose
  int numParticles =100;
  int RNGSeed = 1;
  string outputFileName = "eventGenData.txt";
  for(int iA = 1; iA < argc;iA++)
  {
    if (!strcmp(argv[iA], "-n")) //number of particles
    {
      numParticles = stoi(argv[iA+1]);
    }
    if (!strcmp(argv[iA], "-s")) //RNG seed
    {
      RNGSeed = stoi(argv[iA+1]);
    }
    if (!strcmp(argv[iA], "-o"))
    {
      outputFileName = argv[iA+1];
    }
  }

  extraMath eMath;

  randNum randNumGen;
  randNumGen.setSeed(RNGSeed);
  cout << "Seed for randum number generator = " << randNumGen.m_si<<"\n";
  /// read in simConfig file:
  /// read in target data
  /// read in beam data
  ifstream ifs("../simConfig.dat");
  string line;
  getline(ifs, line);
  string isotope = "26_56_Iron";
  double targetWidth;
  double targetThickness;
  vector<double> beamEner;
  vector<double> beamFlux;
  vector<double> fluxCDF;
  while (line!="#end")
  {
    if (line.substr(0,8)=="#isotope") 
    {
      getline(ifs, line);
      isotope = line;
    }
    else if (line.substr(0,6)=="#width") 
    {
      getline(ifs, line);
      targetWidth = (double)stof(line.substr(0,line.find(" ")));
    }
    else if (line.substr(0,10)=="#thickness") 
    {
      getline(ifs, line);
      targetThickness = (double)stof(line.substr(0,line.find(" ")));
    }
    else if (line.substr(0,5)=="#flux")
    {
      getline(ifs,line);
      if (line=="#hist")
      {
        cout << "Reading custom flux\n";
        getline(ifs,line);
        while (line!="#endflux")
        {
          beamEner.push_back((double)stof(line.substr(0,line.find(","))));
          beamFlux.push_back((double)stof(line.substr(line.find(",")+1)));
          getline(ifs, line);
        }      
      }
      else
      {
        if (line=="#gauss")
        {
          cout << "Generating Gaussian flux\n";
          getline(ifs,line);
          double e_min = stof(line.substr(0,line.find(" ")));
          double e_max = stof(line.substr(line.find(" ")+1));
          getline(ifs,line);
          double mean = stof(line.substr(0,line.find(" ")));
          double std_dev = stof(line.substr(line.find(" ")+1));
          for (int i=0;i<100;i++)
          {
            double ener = (e_max-e_min)/100*i;
            beamEner.push_back(ener);
            beamFlux.push_back(exp(-pow((mean-ener)/(sqrt(2.)*std_dev),2)));
          }
          fluxCDF = eMath.fluxCDF(beamFlux);
        }
        else if (line=="#uniform")
        {
          cout << "Generating Uniform flux\n";
          getline(ifs,line);
          double e_min = stof(line.substr(0,line.find(" ")));
          double e_max = stof(line.substr(line.find(" ")+1));
          getline(ifs,line);
          for (int i=0;i<100;i++)
          {
            beamEner.push_back((e_max-e_min)/100*i + e_min);
            beamFlux.push_back(1.);
            fluxCDF.push_back(1./(e_max-e_min)*((e_max-e_min)/100*i + e_min));
          }
        }
        else if (line=="#expo")
        {
          cout << "Generating Exponential flux\n";
          getline(ifs,line);
          double e_min = stof(line.substr(0,line.find(" ")));
          double e_max = stof(line.substr(line.find(" ")+1));
          getline(ifs,line);
          double decay_const = stof(line);
          for (int i=0;i<100;i++)
          {
            double ener = (e_max-e_min)/100*i;
            beamEner.push_back(ener);
            beamFlux.push_back(exp(decay_const*ener));
          }
        }
        else if (line=="#watt")
        {
          cout << "Generating Watt Fission flux\n";
          getline(ifs,line);
          double e_min = stof(line.substr(0,line.find(" ")));
          double e_max = stof(line.substr(line.find(" ")+1));
          getline(ifs,line);
          double a = stof(line.substr(0,line.find(" ")));
          double b = stof(line.substr(line.find(" ")+1));
          for (int i=0;i<100;i++)
          {
            double ener = (e_max-e_min)/100*i;
            beamEner.push_back(ener);
            beamFlux.push_back(exp(-a*ener)*.5*(exp(sqrt(b*ener))-exp(-sqrt(b*ener))));
          }
          fluxCDF = eMath.fluxCDF(beamFlux);
        }
        else 
        {
          cout << "Invalid flux type. See document for correct input formatting\n";
          cout << "Exiting\n";
          return 0;
        }
      }
    }
    getline(ifs,line);
  }
  cout <<"width = "<< targetWidth << ", thickness = " <<targetThickness << "\n";
  /// convert energy-flux to CDF
  for (int i=0;i<100;i++)
  {
    cout << beamEner[i] << ", " << fluxCDF[i]<<"\n";
  }
  /// initialize output file (separate class?)
  ofstream ofs {outputFileName.c_str()};
  /// initialize data
  readXSec initData;
  cout << "initializing Capture Cross Sections\n";
  initData.readCapture(isotope);
  initData.setRandGenerator(&randNumGen);
  cout << "initializing Total Cross Sections\n";
  initData.initTotal(isotope);
  cout << "initializing Partial Inelastic Cross Sections\n";
  initData.initPartialInel(isotope);
  cout << "initializing RIPL data\n";
  RIPL riplData;
  riplData.readRIPL(.001);
  for (int i=0;i<10;i++) randNumGen.next();
  cout << "Reading elastic angular distributions\n";
  initData.readElasticAngle(isotope);
  int numCap=0;
  int numEl =0;
  int numInel =0;
  int numIso=0;
  map<string, int> iso_counts;
  /// begin main loop
  cout << "beginning main loop\n";
  map<string, double> inelLevels;
  int per = numParticles/10;
  if (numParticles<=10) per = 1;
  for (int iPart =0; iPart < numParticles; iPart++)
  {
    if (iPart%per==0) cout << "working on history " << iPart << "\n";

    // interaction location
    vector<double> init_pos_xy = initPosition(targetWidth, (&randNumGen));
    double init_pos_z = targetThickness*randNumGen.next();
    double beamEnergy =0.;
    int hBin = distance(fluxCDF.begin(),std::lower_bound(fluxCDF.begin(), fluxCDF.end(), randNumGen.next()));
    double rand2 = randNumGen.next();
    // if (hBin==0) beamEnergy = beamEner[0]*rand2;
    // else beamEnergy = beamEner[hBin] + (beamEner[hBin+1]-beamEner[hBin])*rand2;
    beamEnergy = beamEner[hBin+1] + (beamEner[hBin+2]-beamEner[hBin+1])*rand()/2e9;
    // beamEnergy = 19.9*rand()/2e9  + .1;
    ofs << "["<<beamEnergy<<",";
    ofs << "("<<init_pos_xy[0]<<","<<init_pos_xy[1]<<","<<init_pos_z<<"),";
    // if (beamEner[hBin+1]<4.1 && beamEner[hBin+1]>3.9)
    // {
    //   cout << beamEner[hBin+1] << ", " <<(beamEner[hBin+2]-beamEner[hBin+1])*randNumGen.next() << "\n";
    // }
    // cout << "beamEnergy = " << beamEnergy<<"\n";
    /// determine reaction or lack thereof
    string rxnType = initData.pickReaction(beamEnergy*1e6);
    /// sample initial position
    /// based on chosen reaction, generate output particles
    if (rxnType=="cap")
    {
      numCap++;
    }
    else if (rxnType=="el")
    {
      numEl++;
      double mu_LAB = initData.pickElasticAngle(beamEnergy*1e6, randNumGen.next());
      // double mu_LAB = 2.*randNumGen.next()-1;
      double new_energy = eMath.outgoingEnergy(beamEnergy,931,931*56,mu_LAB, 0);
      ofs << "{n,"<<new_energy<<","<<mu_LAB<<"}";
      // cout << beamEnergy << ": e_out = " << e_out << ", ang = " << mu_LAB << "\n";
    }
    else if (rxnType=="inel")
    {
      numInel++;
      string inelLevel_ensdf = initData.pickInelLevel(beamEnergy*1e6);
      // cout << "ensdfLevel = " << inelLevel_ensdf << "\n";
      uint16_t riplLevel = riplData.getLevelNumber(inelLevel_ensdf, (&randNumGen));
      // cout << "riplLevel = " <<riplLevel << "\n";
      double energy = riplData.getLevelEnergy(riplLevel);
      // cout << "energy = " << energy << "\n";

      vector<double> gammas = riplData.doCascade(riplLevel, (&randNumGen));
      double mu_LAB = initData.pickInelAngle(inelLevel_ensdf, beamEnergy*1e6, randNumGen.next());
      double new_energy = eMath.outgoingEnergy(beamEnergy,931,931*56,mu_LAB, energy);
      ofs << "{n,"<<new_energy<<","<<mu_LAB<<"},";
      inelLevels[inelLevel_ensdf] = energy;
      //for (double gamma : gammas)
      for (int iG = 0; iG<gammas.size()-1; iG++)
      {
        ofs <<"{g,"<< gammas[iG] <<","<<2.*randNumGen.next()-1<<"},";
      }
      ofs <<"{g,"<< gammas[gammas.size()-1] <<","<<2.*randNumGen.next()-1<<"}";
    }
    else
    {
      numIso++;
      iso_counts[rxnType]+=1;
      // do decay
      vector<double> gammas = initData.readENSDF(rxnType, randNumGen.next());
      if (gammas.size()>0)
      {
        for (int iG = 0; iG<gammas.size()-1; iG++)
        {
          // cout << "decay gamma = " << gammas[iG] << "\n";
          ofs << "{g," << gammas[iG] << 2.*randNumGen.next()-1<<"},";
        }
        ofs << "{g," << gammas[gammas.size()-1] << 2.*randNumGen.next()-1<<"}";
      }
    }
    ofs << "]\n";
  }

  cout << "numCaptures = " << numCap <<"\n";
  cout << "numElastic = " << numEl <<"\n";
  cout << "numInelastics = " << numInel <<"\n";
  cout << "numIsos = " << numIso << "\n";
  for (auto& pair : inelLevels)
  {
    cout << pair.first << ": " << pair.second << "\n";
  }
  for (auto& pair : iso_counts)
  {
    cout << pair.first << ": " << pair.second << "\n";
  }
  ofs.close();
  return 0;
}