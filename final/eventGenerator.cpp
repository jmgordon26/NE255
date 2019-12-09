// #include "dataTypes.h"
// #include "randNum.h"
#include "readXSec.h"
#include "extraMath.h"

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

using namespace std;

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

  ofstream ofs {outputFileName.c_str()};

  extraMath eMath;

  randNum randNumGen;
  randNumGen.setSeed(RNGSeed);
  cout << "Seed for randum number generator = " << randNumGen.m_si<<"\n";
  /// read in simConfig file:
  /// read in target data
  /// read in beam data
  ifstream ifs("simConfig.dat");
  string line;
  getline(ifs, line);
  string isotope = "26_56_Iron";
  string targetWidth;
  string targetThickness;
  vector<double> beamEner;
  vector<double> beamFlux;
  // while (line!="#end")
  // {
  //   if (line=="#isotope") isotope = line;
  //   else if (line=="#width") targetWidth = line;
  //   else if (line=="#thickness") targetThickness = line;
  //   else if (line=="#flux")
  //   {
  //     getline(ifs,line);
  //     while (line!="#endflux")
  //     {
  //       beamEner.push_back(stof(line.substr(0,line.find(","))));
  //       beamFlux.push_back(stof(line.substr(line.find(",")+1)));
  //       getline(ifs, line);
  //     }
  //   }
  //   getline(ifs,line);
  // }
  /// convert energy-flux to CDF
  vector<double> hHE = {1.046,  1.073,  1.101,  1.131,  1.16 ,  1.192,  1.225,  1.259,
        1.295,  1.333,  1.371,  1.412,  1.455,  1.499,  1.546,  1.594,
        1.646,  1.7  ,  1.756,  1.816,  1.877,  1.943,  2.012,  2.084,
        2.161,  2.243,  2.328,  2.419,  2.515,  2.617,  2.726,  2.841,
        2.963,  3.094,  3.235,  3.384,  3.544,  3.716,  3.9  ,  4.099,
        4.313,  4.545,  4.796,  5.068,  5.365,  5.688,  6.042,  6.429,
        6.855,  7.325,  7.845,  8.423,  9.068,  9.79 , 10.602, 11.521,
       12.565, 13.759, 15.133, 16.726, 18.587, 20.781};
  double hLE = 1.02;
  vector<double> hP = {0.        , 0.00251036, 0.00427961, 0.00580133, 0.00769592,
       0.01082054, 0.01298782, 0.01490305, 0.01740344, 0.02015121,
       0.02279021, 0.02510506, 0.02743855, 0.03079164, 0.03382312,
       0.03776453, 0.04166598, 0.0455366 , 0.04949852, 0.05373237,
       0.05829342, 0.06276764, 0.06703562, 0.07052072, 0.07485722,
       0.08022411, 0.08515162, 0.09093775, 0.09680106, 0.1038973 ,
       0.11221593, 0.12047039, 0.12946396, 0.14089085, 0.1507597 ,
       0.16332836, 0.17605697, 0.18821974, 0.20302229, 0.21768669,
       0.23473363, 0.25438739, 0.27677956, 0.30068547, 0.32798623,
       0.35678972, 0.3905468 , 0.42910989, 0.47566086, 0.5211588 ,
       0.5684784 , 0.62145803, 0.67450177, 0.73019106, 0.78721138,
       0.84650753, 0.89838763, 0.95017052, 0.98117838, 0.99365541,
       0.99891237, 1.};
  /// initialize output file (separate class?)

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
  for (int iPart =0; iPart < numParticles; iPart++)
  {
    if (iPart%1000==0) cout << "working on history " << iPart << "\n";
    /// sample beam energy, location, direction
    // int hBin = (*upper_bound(beamFlux, beamFlux+62, randNumGen.next()));
    // double rand2 = randNumGen.next();
    // if (hBin==1) energy = hLE + (beamEner[0]-hLE)*rand2;
    // else energy = beamEner[hBin] + (beamEner[hBin+1]-beamEner[hBin])*rand2;
    double beamEnergy =0.;
    int hBin = distance(hP.begin(),std::upper_bound(hP.begin(), hP.end(), randNumGen.next()));
    double rand2 = randNumGen.next();
    if (hBin==1) beamEnergy = hLE + (hHE[0]-hLE)*rand2;
    else beamEnergy = hHE[hBin] + (hHE[hBin+1]-hHE[hBin])*rand2;
    ofs << "{"<<beamEnergy<<",";
    // cout << "beamEnergy = " << beamEnergy<<"\n";
    string rxnType = initData.pickReaction(beamEnergy*1e6);
    /// determine reaction or lack thereof

    /// based on chosen reaction, generate output particles
    if (rxnType=="cap")
    {
      numCap++;
    }
    else if (rxnType=="el")
    {
      numEl++;
      double mu_LAB = initData.pickElasticAngle(beamEnergy*1e6, randNumGen.next());
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
      // double new_energy = pow(56./57,2)*(beamEnergy-energy*57./56);


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
    // else if (rxnType=="iso")
    else
    {
      numIso++;
      iso_counts[rxnType]+=1;
      // do decay
    }
    ofs << "}\n";
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