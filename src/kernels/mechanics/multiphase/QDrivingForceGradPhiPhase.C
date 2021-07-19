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

#include "QDrivingForceGradPhiPhase.h"
registerMooseObject("gibbsApp", QDrivingForceGradPhiPhase);

template<>
InputParameters
validParams<QDrivingForceGradPhiPhase>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Mechanical driving traction of phase transformation");
  params.addRequiredParam<std::string>("base_name", "The phase name");
  params.addRequiredParam<MaterialPropertyName>("mob_name", "phase field mobility");
  params.addRequiredParam<MaterialPropertyName>("nd_factor", "Non-dimensional factor");
  return params;
}

QDrivingForceGradPhiPhase::QDrivingForceGradPhiPhase(const InputParameters & parameters)
  :Kernel(parameters),
  _phase_name(getParam<std::string>("base_name")),
  _dfbulk_dgradphi(getMaterialProperty<RealVectorValue>(_phase_name + "_dfbulk_dgradphi")),
  _d2fbulk_dgradphi_dphi(getMaterialProperty<RealVectorValue>(_phase_name + "_d2fbulk_dgradphi_dphi")),
  _L(getMaterialProperty<Real>
             (getParam<MaterialPropertyName>("mob_name"))),
  _nd_factor(getMaterialProperty<Real>
             (getParam<MaterialPropertyName>("nd_factor")))
{
}

Real
QDrivingForceGradPhiPhase::computeQpResidual(){
  return _grad_test[_i][_qp] *_L[_qp] * _dfbulk_dgradphi[_qp] * _nd_factor[_qp];
}

Real
QDrivingForceGradPhiPhase::computeQpJacobian(){
  return _grad_test[_i][_qp] * _L[_qp] * _d2fbulk_dgradphi_dphi[_qp]
                                       * _nd_factor[_qp] * _phi[_j][_qp];
}

Real
QDrivingForceGradPhiPhase::computeQpOffDiagJacobian(unsigned int /*jvar*/){
    return 0;
}
