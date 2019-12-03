#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>
using namespace std;
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

class particle
{
public:
  particle();
  void setLoc(double a_loc) {m_loc = a_loc;}
  double getLoc() {return m_loc;}

  void setDir(double a_dirn) {m_dirn = a_dirn/abs(a_dirn);}
  double getDir() {return m_dirn;}

  bool getStatus() {return m_status;}
  void setStatus(bool a_status) {m_status = a_status;}

  void setRegion(int a_region) {m_region = a_region;}
  int getRegion() {return m_region;}

  void print();
private:
  double m_loc = -10.;
  double m_dirn = -1.;
  bool m_status = true;
  int m_region = -10;
};

class region
{
public:
  region();
  void setScatter(double a_scatter);
  void setAbs(double a_abs);
  void setLims(vector<double> a_lims) {m_lims = a_lims;}
  
  void setRegionID(int a_regionID) {m_regionID = a_regionID;}
  vector<double> getLims() {return m_lims;}
  int checkRegion(particle* a_particle);
  double getNextLoc(double a_rand, particle* a_particle);

  int getColl(double a_rand);
private:
  int m_regionID=0;
  double m_scatter = -1.;
  double m_abs = -1.;

  vector<double> m_lims = {-100,100};

  double m_tot = m_scatter+m_abs;
};