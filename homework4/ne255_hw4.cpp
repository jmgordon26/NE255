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
#include "TGraph.h"
struct randNum
{
  int m_c = 0;
  int m_g = 7;
  uint64_t m_p = pow(2,45);

  int m_s0 = 1;
  uint64_t m_si = m_s0;

  double next();
};

double randNum::next()
{
  m_si = (m_si*m_g + m_c)%m_p;
  return 1.*m_si/m_p;
}

vector<TGraph*> hw4()
{
  randNum r1;

  double x1 = 0;
  double y1 = 0;

  double pi[3];
  double absErr[3];
  double relErr[3];
  double it[3] = {100,1000,10000};

  cout << "using LCPRNG with seed = " << r1.m_s0 << "\n";
  int counter=0;
  for (int numIters : {100,1000,10000})
  {
    int numHits = 0;
    for (int i=0; i<numIters; i++)
    {
      x1 = r1.next();
      y1 = r1.next();
      double r = x1*x1 + y1*y1;
      if (sqrt(r)<=1) numHits++;
    }  
    cout << "numIterations = " << numIters <<": pi = " << 4.*numHits/numIters << "\n";
    pi[counter] = 4.*numHits/numIters;
    absErr[counter] = abs(4.*numHits/numIters-3.14159);
    relErr[counter] = abs((4.*numHits/numIters-3.14159) / 3.14159);
    cout << "\trelErr="<<relErr[counter] << ", absErr="<<absErr[counter]<<"\n";
    counter++;
  }

  vector<TGraph*> outp;
  outp.push_back(new TGraph(3,it,pi));
  outp.push_back(new TGraph(3,it,absErr));
  outp.push_back(new TGraph(3,it,relErr));

  return outp;
}