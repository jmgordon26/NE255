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
// project includes
#include "ne255_hw3.h"


using namespace std;

diamondDiff_1D::diamondDiff_1D()
{
  ;
}
void diamondDiff_1D::buildInitFlux(vector<double> angles)
{
  cout << "NumCells = " << m_numCells << "\n";
  for (double angle : angles)
  {
    // for (int iC =0; iC<m_numCells; iC++)
    // {
      m_centerFlux[angle].push_back(1.);
      
      if (angle>0) 
      {
        m_incomingFlux[angle].push_back(2.);
        m_prevFlux[angle].push_back(2.);
      }
      else 
      {
        m_incomingFlux[angle].push_back(0.);
        m_prevFlux[angle].push_back(0.);
      }
    // }
      for (int iC=1; iC<m_numCells; iC++) {m_prevFlux[angle].push_back(0.);}
  }
}

void diamondDiff_1D::setAlpha(double a_alpha)   {m_alpha = a_alpha;}
void diamondDiff_1D::setSpacing(double a_h)     
{
  m_h = a_h;
  m_numCells = (int)m_range/m_h;
}
void diamondDiff_1D::setSigmaT(double a_sigmaT) {m_sigmaT = a_sigmaT;}
void diamondDiff_1D::setSigmaS(double a_sigmaS) {m_sigmaS = a_sigmaS;}
void diamondDiff_1D::setSource(double a_q)      {m_q = a_q;}

void diamondDiff_1D::performForwardSweep(int a_index)
{

  for (auto& pair : m_centerFlux)
  {

    double incomingFlux = *(m_incomingFlux[pair.first].end()-1);
    if (pair.first>0) 
    {
      double fluxTilde =0.;
      for (auto& pair : m_prevFlux)
      {
        fluxTilde+=.5*pair.second[a_index];
        // if (pair.first>0) fluxTilde+=.5*pair.second[a_index];
        // else fluxTilde += .5*pair.second[pair.second.size()-1-a_index];
      }
      double source = 2*3.14159*fluxTilde*m_sigmaS + m_q;
      cout << "source = " << source<<"\n";

      double alpha_denom = 1./(1+m_alpha);
      double cellFlux = source + 2.*alpha_denom*abs(pair.first)/m_h*incomingFlux/(m_sigmaT + 2.*alpha_denom*abs(pair.first)/m_h);
      if (a_index==0) cellFlux=2.0;
      pair.second.push_back(cellFlux);
      // cout << cellFlux <<", ";
      double halfFlux = 2.*alpha_denom*cellFlux - (1-m_alpha)*alpha_denom*incomingFlux;
      if (a_index==m_numCells-1) 
      {
        cout << halfFlux<<"\n";
        cout << 2*alpha_denom*cellFlux - (1+m_alpha)/(1 - m_alpha)*halfFlux << ", ";
        cout << *(m_incomingFlux[pair.first].end()-1) << ",,  ";
      }
      m_incomingFlux[pair.first].push_back(halfFlux);
      cout << pair.first << ", " << alpha_denom << ", " << *(m_incomingFlux[pair.first].end()-1) << "\n";

    }
    else 
    {
      double fluxTilde =0.;
      for (auto& pair : m_prevFlux)
      {
        fluxTilde += .5*pair.second[pair.second.size()-1-a_index];
      }
      double source = 2*3.14159*fluxTilde*m_sigmaS + m_q;
      cout << "source = " << source<<"\n";
      double incomingFlux = m_incomingFlux[pair.first][0];
      double alpha_denom = 1./(1-m_alpha);

      double cellFlux = source + 2.*alpha_denom*abs(pair.first)/m_h*incomingFlux/(m_sigmaT + 2.*alpha_denom*abs(pair.first)/m_h);
      // m_centerFlux[a_angle].insert(m_centerFlux[a_angle].begin(), cellFlux);
      pair.second.insert(pair.second.begin(), cellFlux);

      double halfFlux = 2.*alpha_denom*cellFlux - (1+m_alpha)*alpha_denom*incomingFlux;
      m_incomingFlux[pair.first].insert(m_incomingFlux[pair.first].begin(), halfFlux);
      cout << pair.first << ", " << alpha_denom << ", " << m_incomingFlux[pair.first][0] <<", ";
      cout << incomingFlux <<", "<< 2.*alpha_denom*abs(pair.first)/m_h*incomingFlux/(m_sigmaT + 2.*alpha_denom*abs(pair.first)/m_h)<<"\n ";

    }
    
  }
}

void diamondDiff_1D::performBackwardSweep(int a_index)
{
  double fluxTilde =0.;
  for (auto& pair : m_prevFlux)
  {
    fluxTilde+=pair.second[a_index];
  }
  double source = 2*3.14159*fluxTilde*m_sigmaS + m_q;
  for (auto& pair : m_centerFlux)
  {
    double incomingFlux = m_incomingFlux[pair.first][0];
    double alpha_denom = 1./(1-m_alpha);

    double cellFlux = m_q + 2.*alpha_denom*abs(pair.first)/m_h*incomingFlux/(m_sigmaT + 2.*alpha_denom*abs(pair.first)/m_h);
    // m_centerFlux[a_angle].insert(m_centerFlux[a_angle].begin(), cellFlux);
    pair.second.insert(pair.second.begin(), cellFlux);

    double halfFlux = 2.*alpha_denom*cellFlux - (1+m_alpha)*alpha_denom*incomingFlux;
    m_incomingFlux[pair.first].insert(m_incomingFlux[pair.first].begin(), halfFlux);
      
  }
  
}

double diamondDiff_1D::calcError()
{
  vector<double> errors;
  for (auto& pair : m_centerFlux)
  {
    for (int iFl =0; iFl<pair.second.size(); iFl++)
    {
      errors.push_back(abs(pair.second[iFl] - m_prevFlux[pair.first][iFl]));
    }
  }


  sort(errors.begin(),errors.end());

  m_prevFlux = m_centerFlux;
  for (auto& pair : m_centerFlux)
  {
    pair.second.clear();
    m_incomingFlux[pair.first].clear();
    if (pair.first>0) m_incomingFlux[pair.first].push_back(2.);
    else m_incomingFlux[pair.first].push_back(0.);
  }
  return *(errors.end()-1);
}
