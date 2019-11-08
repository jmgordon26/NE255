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

using namespace std;

class diamondDiff_1D {
public:
  /// default constructor
  diamondDiff_1D();
private:

  // flux at each angle mu
  map<double, vector<double>> m_centerFlux;
  map<double, vector<double>> m_prevFlux;
  map<double, vector<double>> m_incomingFlux;

  double m_alpha = 0.;
  double m_sigmaT = 1.;
  double m_sigmaS = 0.;
  double m_h = .08;
  double m_q = 0.;
  double m_angWeight = .5;
  double m_range = 2.;
public:
  // setters
  void setAlpha(double a_alpha); 
  void setSpacing(double a_h);
  void setSigmaT(double a_sigmaT);
  void setSigmaS(double a_sigmaS);
  void setSource(double a_q);   

  // initial flux to zeros
  void buildInitFlux(vector<double> angles);

  // perform sweep
  void performForwardSweep(int a_index);
  void performBackwardSweep(int a_index);

  // calculate error of current iteration
  // currently does uniform norm
  double calcError();

  int m_numCells = (int)m_range/m_h;
  
};
