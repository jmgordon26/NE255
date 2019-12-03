#include "hw5.h"
#include<iostream>
particle::particle()
{;}
void particle::print()
{
  cout << "Location = " << m_loc;
  cout << "\nStatus = " << m_status;
  cout << "\nDir = " << m_dirn;
  cout << "\nRegion = " << m_region << "\n";;
}
region::region()
{;}
void region::setScatter(double a_scatter) 
{ 
  m_scatter=a_scatter; 
  m_tot = m_scatter+m_abs;
}

void region::setAbs(double a_abs) 
{
  m_abs = a_abs; 
  m_tot = m_scatter+m_abs;
}

double region::getNextLoc(double a_rand, particle* a_particle)
{
  double new_loc = a_particle->getLoc() + a_particle->getDir()*-log(a_rand)/(m_tot);
  return new_loc;
}
int region::checkRegion(particle* a_particle)
{
  if (a_particle->getLoc()<m_lims[1] && a_particle->getLoc()>m_lims[0]) return m_regionID;
  else if (a_particle->getLoc()>=m_lims[1]) return m_regionID+1;
  return m_regionID-1;
}
int region::getColl(double a_rand)
{
  if (a_rand<= m_abs/m_tot) return 0;
  return 1;
}
