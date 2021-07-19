//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html 

//* This kernel implements the chemical potential
//* The equation this kernel mu - df/dc = 0
//* mu is the variable that this kernel operates on

#include "DrivingForceWBM.h"

registerMooseObject("gibbsApp", DrivingForceWBM);

template<>
InputParameters
validParams<DrivingForceWBM>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn:dh*[f_beta -f_alpha] = 0");
  params.addRequiredCoupledVar("c", "composition field of the solute component");
  params.addRequiredParam<MaterialPropertyName>("mob_name", "phase field mobility");
  params.addRequiredParam<MaterialPropertyName>("nd_factor", "Dimensionless factor = RT/(m*Vm)");
  return params;
}

DrivingForceWBM::DrivingForceWBM(const InputParameters & parameters)
  :Kernel(parameters),
  _comp(coupledValue("c")),
  _comp_var(coupled("c")),
  _dh(getMaterialProperty<Real>("dh")),
  _d2h(getMaterialProperty<Real>("d2h")),
  _free_energy_alpha(getMaterialProperty<Real>("f_alpha")),
  _dfalpha_dc(getMaterialProperty<Real>("B_diff_pot_alpha")),
  _free_energy_beta(getMaterialProperty<Real>("f_beta")),
  _dfbeta_dc(getMaterialProperty<Real>("B_diff_pot_beta")),
  _L(getMaterialProperty<Real>("mob_name")),
  _nd_factor(getMaterialProperty<Real>("nd_factor"))
{
}

Real
DrivingForceWBM::computeQpResidual()
{
  //Since the model is WBM
  //We have two diffusion potentials for each phase
  //and only one composition variable
  //Residual: N_{i} * h_{\phi} (f_beta - f_alpha) = 0 ;
  return _test[_i][_qp] *_L[_qp] * _nd_factor[_qp]* _dh[_qp] * (_free_energy_beta[_qp] - _free_energy_alpha[_qp]);

}

Real
DrivingForceWBM::computeQpJacobian()
{
  // This will be non-zero since the 
  // free energies of each phase is dependent on composition
  return  _test[_i][_qp] * _L[_qp] *_nd_factor[_qp]* _d2h[_qp] * (_free_energy_beta[_qp] - _free_energy_alpha[_qp]) *_phi[_j][_qp];
}

Real
DrivingForceWBM::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _comp_var)
 {
    return _test[_i][_qp] *_L[_qp] * _nd_factor[_qp]* _dh[_qp] *(_dfbeta_dc[_qp] - _dfalpha_dc[_qp]) * _phi[_j][_qp];  
 }
 else
 {
    return 0.0;
 }
}
