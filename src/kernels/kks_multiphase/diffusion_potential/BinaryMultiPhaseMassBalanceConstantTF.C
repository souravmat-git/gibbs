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

#include "BinaryMultiPhaseMassBalanceConstantTF.h"

registerMooseObject("gibbsApp", BinaryMultiPhaseMassBalanceConstantTF);

template <>
InputParameters
validParams<BinaryMultiPhaseMassBalanceConstantTF>()
{
  InputParameters params = validParams<BinaryMultiPhaseMassBalance>();
  return params; 
}

BinaryMultiPhaseMassBalanceConstantTF::BinaryMultiPhaseMassBalanceConstantTF(const InputParameters & parameters)
  : BinaryMultiPhaseMassBalance(parameters)
{
}

Real
BinaryMultiPhaseMassBalanceConstantTF::computeQpResidual()
{   
  //Variable on which the kernel operates: X_B
  return (BinaryMultiPhaseMassBalance::computeQpResidual());
}

Real
BinaryMultiPhaseMassBalanceConstantTF::computeQpJacobian()
{
  return (_grad_test[_i][_qp] * ((BinaryMultiPhaseMassBalance::L_BB_interp() * 
                                  BinaryMultiPhaseMassBalance::thermodynamic_factor() * _grad_phi[_j][_qp])
                                +(BinaryMultiPhaseMassBalance::dL_BB_muB_interp()* 
                                  BinaryMultiPhaseMassBalance::thermodynamic_factor()*_grad_u[_qp]* _phi[_j][_qp])) );
}

Real
BinaryMultiPhaseMassBalanceConstantTF::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _B_diff_pot_var)
  { 
    return (BinaryMultiPhaseMassBalance::computeQpOffDiagJacobian(_B_diff_pot_var));
  }
 else if (jvar == _phase_alpha_var)
 {
    return (BinaryMultiPhaseMassBalance::computeQpOffDiagJacobian(_phase_alpha_var));
 }
 else if (jvar == _phase_beta_var)
 {
   return (BinaryMultiPhaseMassBalance::computeQpOffDiagJacobian(_phase_beta_var));
 }
 else
 {
    return 0;
 }
}
