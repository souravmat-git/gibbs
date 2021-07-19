//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* This was written by S.Chatterjee

#include "ContinuityEquation.h"

registerMooseObject("gibbsApp", ContinuityEquation);

template <>
InputParameters
validParams<ContinuityEquation>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("This kernel implements the continuity equation"
                              "Eqn:  \nabla * (M(\nabla mu) = 0");
  params.addRequiredCoupledVar("eta", "To distinguish phases");
  params.addRequiredCoupledVar("xB",  "Mole fraction of component B");
  return params;
}

ContinuityEquation::ContinuityEquation(const InputParameters & parameters)
  : Kernel(parameters),
  _eta(coupledValue("eta")),
  _eta_var(coupled("eta")),
  _xB(coupledGradient("xB")),
  _xB_var(coupled("xB")),
  //thermodynamic_factor
  _B_tf_alpha(getMaterialProperty<Real>("B_therm_factor_alpha")),
  _B_tf_beta(getMaterialProperty<Real>("B_therm_factor_beta")),
  _h(getMaterialProperty<Real>("h")),
  _dh(getMaterialProperty<Real>("dh")),
  //Kinetic materials
  _L_BB_alpha(getMaterialProperty<Real>("L_BB_alpha")),
  _L_BB_beta(getMaterialProperty<Real>("L_BB_beta")),
  _dL_BB_xB_alpha(getMaterialProperty<Real>("dL_BB_xB_alpha")),
  _dL_BB_xB_beta(getMaterialProperty<Real>("dL_BB_xB_beta"))
{
}

Real
ContinuityEquation::thermodynamic_factor() const{
  return (_h[_qp] * _B_tf_alpha[_qp] + (1.0 -_h[_qp])* _B_tf_alpha[_qp]);
}

Real
ContinuityEquation::L_BB_interp() const{
  return (_h[_qp] * _L_BB_beta[_qp] + (1-_h[_qp])* _L_BB_alpha[_qp]);
}

Real
ContinuityEquation::dL_BB_xB_interp() const{
  return (_h[_qp] * _dL_BB_xB_beta[_qp] + (1-_h[_qp])* _dL_BB_xB_alpha[_qp]);
}

Real
ContinuityEquation::computeQpResidual()
{
  return (_grad_test[_i][_qp] * ContinuityEquation::L_BB_interp() * _grad_u[_qp]);
}

Real
ContinuityEquation::computeQpJacobian()
{
  //Derivative wrt B_diff_potential
  return (_grad_test[_i][_qp] * ContinuityEquation::L_BB_interp() * _grad_phi[_j][_qp]);
}

Real
ContinuityEquation::computeQpOffDiagJacobian(unsigned int jvar)
{ 

  if (jvar == _eta_var)
  {
    return (_grad_test[_i][_qp]* _dh[_qp]*(_L_BB_beta[_qp] -_L_BB_alpha[_qp])*_grad_u[_qp]*_phi[_j][_qp]);
  }
  else if (jvar == _xB_var)
  {
    return (_grad_test[_i][_qp]* ContinuityEquation::L_BB_interp()*ContinuityEquation::thermodynamic_factor()*_grad_phi[_j][_qp]
          + _grad_test[_i][_qp]* ContinuityEquation::dL_BB_xB_interp()*_grad_u[_qp]*_phi[_j][_qp]);
  }
  else
  {
    return 0;
  }    
}
