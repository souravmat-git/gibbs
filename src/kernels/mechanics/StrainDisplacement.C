//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "StrainDisplacement.h"
registerMooseObject("gibbsApp", StrainDisplacement);

template<>
InputParameters
validParams<StrainDisplacement>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Strain displacement relation");
  params.addRequiredCoupledVar("ex", "Strain in the x-direction");
  return params;
}

StrainDisplacement::StrainDisplacement(const InputParameters & parameters)
 : Kernel(parameters),
  _ex(coupledValue("ex")),
  _ex_var(coupled("ex"))
{
}

Real
StrainDisplacement::computeQpResidual(){  
//The non-linear variable that this kernel acts on is the displacement variable
  return _test[_i][_qp] * _ex[_qp]  
       + _grad_test[_i][_qp](0)*_u[_qp];
}

Real
StrainDisplacement::computeQpJacobian(){
  return _grad_test[_i][_qp](0) * _phi[_j][_qp];
}

Real
StrainDisplacement::computeQpOffDiagJacobian(unsigned int jvar){
 if (jvar == _ex_var)
  return _test[_i][_qp] * _phi[_j][_qp];
 else 
  return 0;
}
