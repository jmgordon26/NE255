#include "extraMath.h"
#include <stdio.h>
#include <iostream>
double extraMath::outgoingEnergy(double a_incEnergy, double a_mN, double a_mT, double a_labAngle, double a_extraE)
// double extraMath::outgoingEnergy(double a_incEnergy,  double a_inMass, double a_tarMass, double a_labAngle)
{
  double mn = a_mN*1e-13/pow(3e8,2);
  double mt = a_mT*1e-13/pow(3e8,2);
  double gamma = a_incEnergy/a_mN + 1;
  double vel = 3e8*sqrt(1-1./pow(gamma,2));
  double p0 = gamma*mn*vel;
  // double e0 = gamma*mn*pow(3e8,2);

  // double a = pow(mn+mt,2)*pow(3e8,4) + e0*e0 + 2*e0*(mn + mt)*pow(3e8,2) - p0*p0*pow(3e8,2)*a_labAngle*a_labAngle;
  // double b = -2.*e0*(mt*mt - mn*mn)*pow(3e8,4) - 2*e0*e0*(mt - mn)*pow(3e8,2);
  // double c = -1.*(p0*p0*pow(3e8,2)*mn*mn*a_labAngle*a_labAngle*pow(3e8,4) - e0*e0*pow(mt - mn,2)*pow(3e8,4));

  double p0c = gamma*a_mN*vel/3e8;
  double e0 = gamma*a_mN;
  // double a = pow(a_mT+a_mN+e0,2)-p0c*p0c*a_labAngle*a_labAngle;
  // double b = -1.*(e0*(a_mT+a_mN+e0)*(a_mT-a_mN));
  // double c = e0*e0*pow(a_mT-a_mN,2) + p0c*p0c*a_mN*a_mN*a_labAngle*a_labAngle;
  // double a = pow(a_mN+a_mT,2) + e0*e0 + 2*e0*(a_mN + a_mT) - p0c*p0c*a_labAngle*a_labAngle;
  // double b = -2.*e0*(a_mT*a_mT - a_mN*a_mN) - 2*e0*e0*(a_mT - a_mN);
  // double c = (p0c*p0c*a_mN*a_mN*a_labAngle*a_labAngle + e0*e0*pow(a_mT - a_mN,2));

  // double e_out_pos = (-b + sqrt(b*b - 4.*a*c))/(2*a);
  // double e_out_neg = (-b - sqrt(b*b - 4.*a*c))/(2*a);
  // std::cout << "a = "<<a<<", b = " <<b<<", c = "<<c<<"\n";
  // std::cout << p0c << " vs " << e0 << ": " << e_out_pos << " vs " << e_out_neg << "\n";
  // if (e_out_pos>0) return e_out_pos-a_mN;
  // else return e_out_neg-a_mN;

  // double e_in = a_incEnergy + a_inMass;
  // double alpha = pow(a_inMass+a_tarMass,2) + e_in*e_in + 2*e_in*(a_inMass+a_tarMass) - pow(a_incEnergy,2)*pow(a_labAngle,2);
  // double beta = -2*e_in*(a_inMass*a_inMass - a_tarMass*a_tarMass) -2*e_in*e_in*(a_tarMass - a_inMass);
  // double gamma = -1.*(a_incEnergy*a_incEnergy*a_inMass*a_inMass*a_labAngle*a_labAngle + e_in*e_in*pow(a_tarMass- a_inMass,2));
  // double e_out_pos = (-beta + sqrt(beta*beta - 4.*alpha*gamma))/(2*alpha);
  // double e_out_neg = (-beta - sqrt(beta*beta - 4.*alpha*gamma))/(2*alpha);
  // if (e_out_pos>0) return e_out_pos;
  // else return e_out_neg;

  p0 = sqrt(2*1.67e-27*a_incEnergy*1.6e-13);
  double a = 1+a_mN/a_mT;
  double b = 2*a_mN/a_mT*p0*a_labAngle;
  double c = p0*p0*(a_mN/a_mT-1) + 2.*1.67e-27*a_extraE*1.6e-13;

  double e_out_pos = pow((-b + sqrt(b*b - 4.*a*c))/(2*a),2)/(2*1.67e-27)/1.6e-13;
  double e_out_neg = pow((-b - sqrt(b*b - 4.*a*c))/(2*a),2)/(2*1.67e-27)/1.6e-13;
  // std::cout << "a = "<<a<<", b = " <<b<<", c = "<<c<<"\n";
  // std::cout << p0 << " vs " << e0 << ": " << e_out_pos << " vs " << e_out_neg << "\n";
  if (e_out_pos>0) return e_out_pos;
  else return e_out_neg;

} 

std::vector<double> extraMath::outgoingAngle(std::vector<double> a_initAngle, double a_muLAB)
{
  return {a_initAngle[0]+a_muLAB, a_initAngle[1], a_initAngle[2]};
}