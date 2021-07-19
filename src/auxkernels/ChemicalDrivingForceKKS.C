#include "ChemicalDrivingForceKKS.h"

registerMooseObject("gibbsApp",ChemicalDrivingForceKKS);

template<>
InputParameters
validParams<ChemicalDrivingForceKKS>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("diff_pot", "Phase composition in alpha phase");
  
  return params;
}

ChemicalDrivingForceKKS::ChemicalDrivingForceKKS(const InputParameters & parameters)
  :AuxKernel(parameters),
  _diff_pot(coupledValue("diff_pot")),
  _dh(getMaterialProperty<Real>("dh")),
  _omega_alpha(getMaterialProperty<Real>("omega_alpha")),
  _omega_beta(getMaterialProperty<Real>("omega_beta"))
{
}
 
Real
ChemicalDrivingForceKKS::computeValue()
{
  
  return (_omega_beta[_qp] - _omega_alpha[_qp])*_dh[_qp];
}
