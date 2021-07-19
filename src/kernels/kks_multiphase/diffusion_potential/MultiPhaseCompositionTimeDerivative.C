//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MultiPhaseCompositionTimeDerivative.h"

registerMooseObject("gibbsApp",MultiPhaseCompositionTimeDerivative);

template <>
InputParameters
validParams<MultiPhaseCompositionTimeDerivative>()
{
  InputParameters params = validParams<CoupledTimeDerivative>();
  params.addClassDescription("A modified coupled time derivative Kernel that multiplies the time "
                             "derivative of a coupled variable by a generalized susceptibility");
  params.addRequiredCoupledVar("phase_2", "Phase_2, Note: (Phase2 - Phase1)");
  params.addCoupledVar("phase_3", 0.0, "Phase_3, Note: (Phase3 - Phase1)");
  params.addRequiredParam<MaterialPropertyName>("xB_1", "Grand-potential phase 1");
  params.addRequiredParam<MaterialPropertyName>("xB_2", "Grand-potential phase 2");
  params.addRequiredParam<MaterialPropertyName>("inv_B_tf_1", "Thermodynamic factor 1");
  params.addRequiredParam<MaterialPropertyName>("inv_B_tf_2", "Thermodynamic factor 2");
  params.addRequiredParam<MaterialPropertyName>("dh", "Derivative of phase2 wr.t 1");
  params.addRequiredParam<MaterialPropertyName>("d2h","Derivative of dh w.r.t 1");
  params.addRequiredParam<MaterialPropertyName>("d2h_2", "Derivative of dh w.r.t 2");
  params.addParam<MaterialPropertyName>("d2h_3", 0.0,"Derivative of dh w.r.t 3");
  return params;
}

MultiPhaseCompositionTimeDerivative::MultiPhaseCompositionTimeDerivative(const InputParameters & parameters)
  :CoupledTimeDerivative(parameters),
  _phase_2(coupledValue("phase_2")),
  _phase_2_var(coupled("phase_2")),
  //Extend to multiple phases
  _phase_3(coupledValue("phase_3")),
  _phase_3_var(coupled("phase_3")),
  //Get all the declared names for material property
  _xB_1_name(getParam<MaterialPropertyName>("xB_1")),
  _xB_2_name(getParam<MaterialPropertyName>("xB_2")),
  _inv_B_tf_1_name(getParam<MaterialPropertyName>("inv_B_tf_1")),
  _inv_B_tf_2_name(getParam<MaterialPropertyName>("inv_B_tf_2")),
  //Get all the declared names for interpolation function
  _dh_name(getParam<MaterialPropertyName>("dh")),
  _d2h_name(getParam<MaterialPropertyName>("d2h")),
  _d2h_2_name(getParam<MaterialPropertyName>("d2h_2")),
  _d2h_3_name(getParam<MaterialPropertyName>("d2h_3")),
  //Material specific properties
  _xB_1(getMaterialProperty<Real>(_xB_1_name)),
  _xB_2(getMaterialProperty<Real>(_xB_2_name)), 
  _inv_B_tf_1(getMaterialProperty<Real>(_inv_B_tf_1_name)),
  _inv_B_tf_2(getMaterialProperty<Real>(_inv_B_tf_2_name)),
  //interpolation function
  _dh(getMaterialProperty<Real>(_dh_name)),
  _d2h(getMaterialProperty<Real>(_d2h_name)),
  _d2h_2(getMaterialProperty<Real>(_d2h_2_name)),
  _d2h_3(getMaterialProperty<Real>(_d2h_3_name))
{
}

Real
MultiPhaseCompositionTimeDerivative::computeQpResidual()
{
  return (CoupledTimeDerivative::computeQpResidual() * (_xB_2[_qp] - _xB_1[_qp])*_dh[_qp]);
}

Real
MultiPhaseCompositionTimeDerivative::computeQpJacobian()
{
  return CoupledTimeDerivative::computeQpResidual() * (_inv_B_tf_2[_qp] - _inv_B_tf_1[_qp]) * _dh[_qp] * _phi[_j][_qp];
}

Real
MultiPhaseCompositionTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _v_var)
  {
    return (CoupledTimeDerivative::computeQpOffDiagJacobian(jvar) * (_xB_2[_qp] - _xB_1[_qp])*_dh[_qp] +
           CoupledTimeDerivative::computeQpResidual() * (_xB_2[_qp] - _xB_1[_qp])*_d2h[_qp] * _phi[_j][_qp]);
  }
  else if (jvar == _phase_2_var)
  {
    return (CoupledTimeDerivative::computeQpResidual() * (_xB_2[_qp] - _xB_1[_qp]) * _d2h_2[_qp] * _phi[_j][_qp]);
  }
  else if (jvar == _phase_3_var)
  {
    return (CoupledTimeDerivative::computeQpResidual() * (_xB_2[_qp] - _xB_1[_qp]) * _d2h_3[_qp] * _phi[_j][_qp]);
  }
  else
  return 0;
}
