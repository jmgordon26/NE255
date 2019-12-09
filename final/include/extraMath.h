#include <cmath>
#include <vector>
class extraMath
{
public:
  double outgoingEnergy(double a_incEnergy, double a_inMass, double a_tarMass, double a_labAngle, double a_extraE);
  std::vector<double> outgoingAngle(std::vector<double> a_initAngle, double a_muLAB);//, randNumGen*  a_randGen);
  std::vector<double> fluxCDF(std::vector<double> a_fluxPDF);
  // std::vector<double> initPosition(double a_rad, randNum* a_randGen);
};