#include "TotalEnergy.h"

registerMooseObject("gibbsApp", TotalEnergy);

template<>
InputParameters
validParams<TotalEnergy>()
{
  InputParameters params = validParams<InterfaceEnergy>();
  params.addRequiredParam<MaterialPropertyName>("f_alpha", "Free energy of the alpha phase");
  params.addRequiredParam<MaterialPropertyName>("f_beta", "Free energyof the beta phase");
  //params.addRequiredParam<MaterialPropertyName>("h_alpha","Interpolation of alpha phase");
  params.addRequiredParam<MaterialPropertyName>("h_beta", "Interpolation of beta phase");
  return params;
}

TotalEnergy::TotalEnergy(const InputParameters & parameters)
  :InterfaceEnergy(parameters),
  _f_alpha(getMaterialProperty<Real>("f_alpha")),
  _f_beta(getMaterialProperty<Real>("f_beta")),
  //_h_alpha(getMaterialProperty<Real>("h_alpha")),
  _h_beta(getMaterialProperty<Real>("h_beta"))
{
}
 
Real
TotalEnergy::computeValue(){   
  Real _h_alpha = (1.0 - _h_beta[_qp]);                         
  return (InterfaceEnergy::computeValue() + 
               _h_beta[_qp] * _f_beta[_qp] + _h_alpha * _f_alpha[_qp]);
}
