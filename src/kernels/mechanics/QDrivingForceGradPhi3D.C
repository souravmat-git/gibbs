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

#include "QDrivingForceGradPhi3D.h"
registerMooseObject("gibbsApp", QDrivingForceGradPhi3D);

template<>
InputParameters
validParams<QDrivingForceGradPhi3D>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Mechanical driving traction of phase transformation");
  //params.addRequiredCoupledVar("ux", "Displacecment in the x-direction");
  //params.addRequiredParam<MaterialPropertyName>("strain_jump_name",
  //                                          "Name of strain jump material");
  params.addRequiredParam<MaterialPropertyName>("mob_name", "phase field mobility");
  params.addRequiredParam<MaterialPropertyName>("nd_factor", "Non-dimensional factor");
  return params;
}

QDrivingForceGradPhi3D::QDrivingForceGradPhi3D(const InputParameters & parameters)
  :Kernel(parameters),
  _df_dgradphi(getMaterialProperty<RealVectorValue>("df_dgradphi")),
  _d2f_dgradphi_dphi(getMaterialProperty<RealVectorValue>("d2f_dgradphi_dphi")),
  _L(getMaterialProperty<Real>
             (getParam<MaterialPropertyName>("mob_name"))),
  _nd_factor(getMaterialProperty<Real>
             (getParam<MaterialPropertyName>("nd_factor")))
{
}

Real
QDrivingForceGradPhi3D::computeQpResidual(){
  return _grad_test[_i][_qp] *_L[_qp] * _df_dgradphi[_qp] * _nd_factor[_qp];
}

Real
QDrivingForceGradPhi3D::computeQpJacobian(){
  return _grad_test[_i][_qp] * _L[_qp] * _d2f_dgradphi_dphi[_qp] *
                                         _nd_factor[_qp] * _phi[_j][_qp];
}

Real
QDrivingForceGradPhi3D::computeQpOffDiagJacobian(unsigned int /*jvar*/){
    return 0;
}
