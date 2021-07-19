#include "ChemicalDrivingForceWBM.h"

registerMooseObject("gibbsApp",ChemicalDrivingForceWBM);

template<>
InputParameters
validParams<ChemicalDrivingForceWBM>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("c", "Composition of the solute component");  
  return params;
}

ChemicalDrivingForceWBM::ChemicalDrivingForceWBM(const InputParameters & parameters)
  :AuxKernel(parameters),
  _comp(coupledValue("c")),
  _dh(getMaterialProperty<Real>("dh")),
  _free_energy_alpha(getMaterialProperty<Real>("f_alpha")),
  _free_energy_beta(getMaterialProperty<Real>("f_beta")), 
  _dfalpha_dc(getMaterialProperty<Real>("B_diff_pot_alpha")),
  _dfbeta_dc(getMaterialProperty<Real>("B_diff_pot_beta"))
{
}
 
Real
ChemicalDrivingForceWBM::computeValue()
{
 
  const Real omega_beta = _free_energy_beta[_qp] - _dfbeta_dc[_qp] * _comp[_qp];
  const Real omega_alpha = _free_energy_alpha[_qp] - _dfalpha_dc[_qp] * _comp[_qp];
  
  return (omega_beta - omega_alpha)* _dh[_qp];
}
