#include<cmath>
#include <cstdint>
using namespace std;

class randNum
{
public:
  randNum() {;}
  int m_c = 0;
  int m_g = 7;
  uint64_t m_p = pow(2,45);

  int m_s0 = 1;
  uint64_t m_si = m_s0;

  double next()
  {
    m_si = (m_si*m_g + m_c)%m_p;
    return 1.*m_si/m_p;
  }
  void setSeed(int a_seed) {m_s0=a_seed; m_si=m_s0;}
};

