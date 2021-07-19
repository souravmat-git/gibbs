//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html 

//* This kernel implements the driving force 
//* which is the difference in elastic strain energy between the two phases

#include "QDrivingForce3DOffDiag.h"
registerMooseObject("gibbsApp", QDrivingForce3DOffDiag);

template<>
InputParameters
validParams<QDrivingForce3DOffDiag>()
{
  InputParameters params = validParams<AllenCahnElasticEnergyOffDiag>();
  params.addRequiredParam<MaterialPropertyName>("nd_factor", "Non-dimensional factor");
  return params;
}

QDrivingForce3DOffDiag::QDrivingForce3DOffDiag(const InputParameters & parameters)
  :AllenCahnElasticEnergyOffDiag(parameters),
  _ndisp(coupledComponents("displacements")),
  _nd_factor(getMaterialProperty<Real>
             (getParam<MaterialPropertyName>("nd_factor")))
{
}

Real
QDrivingForce3DOffDiag::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_ndisp != 0){
    return AllenCahnElasticEnergyOffDiag::computeQpOffDiagJacobian(jvar)*_nd_factor[_qp];
  }
  else
    return 0.0;
}
