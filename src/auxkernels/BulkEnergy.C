#include "BulkEnergy.h"

registerMooseObject("gibbsApp",BulkEnergy);

template<>
InputParameters
validParams<BulkEnergy>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<MaterialPropertyName>("f_alpha","Free energy of the alpha phase");
  params.addRequiredParam<MaterialPropertyName>("f_beta", "Free energy of the beta phase");
  //params.addParam<MaterialPropertyName>("h_alpha","Interpolation of alpha phase");
  params.addRequiredParam<MaterialPropertyName>("h_beta", "Interpolation of beta phase");
  return params;
}

BulkEnergy::BulkEnergy(const InputParameters & parameters)
  :AuxKernel(parameters),
  //_h_alpha(getMaterialProperty<Real>("h_alpha")),
  _h_beta(getMaterialProperty<Real>("h_beta")),
  _f_alpha(getMaterialProperty<Real>("f_alpha")),
  _f_beta(getMaterialProperty<Real>("f_beta"))
{
}
 
Real
BulkEnergy::computeValue(){

  Real _h_alpha = (1 - _h_beta[_qp]);
  return (_h_beta[_qp] * _f_beta[_qp] + _h_alpha * _f_alpha[_qp]);
}
