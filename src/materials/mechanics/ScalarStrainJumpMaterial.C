//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//

#include "ScalarStrainJumpMaterial.h"
registerMooseObject("gibbsApp", ScalarStrainJumpMaterial);

//template <>
InputParameters
ScalarStrainJumpMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredCoupledVar("ux", "Displacement variable");
  params.addRequiredParam<MaterialPropertyName>("mat_const_alpha", "Material constant for alpha phase");
  params.addRequiredParam<MaterialPropertyName>("mat_const_beta", "Material constant for alpha phase");
  params.addRequiredParam<MaterialPropertyName>("eT_alpha", "Transformation strain for alpha phase");
  params.addRequiredParam<MaterialPropertyName>("eT_beta", "Transformation strain for beta phase");
  params.addRequiredParam<MaterialPropertyName>("h", "Interpolation Function");
  params.addRequiredParam<MaterialPropertyName>("a","Jump in strain");
  params.addRequiredParam<MaterialPropertyName>("da_de", "First derivative of jump in strain wrt overall strain");
  params.addRequiredParam<MaterialPropertyName>("da_dh", "First derivative of jump in strain wrt h");
  return params;
}

ScalarStrainJumpMaterial::ScalarStrainJumpMaterial(const InputParameters & parameters)
  : Material(parameters),
   _grad_ux(coupledGradient("ux")),
   _mat_const_alpha(getMaterialProperty<Real>(getParam<MaterialPropertyName>("mat_const_alpha"))),
   _mat_const_beta(getMaterialProperty<Real>(getParam<MaterialPropertyName>("mat_const_beta"))),
   _eT_alpha(getMaterialProperty<Real>(getParam<MaterialPropertyName>("eT_alpha"))),
   _eT_beta(getMaterialProperty<Real>(getParam<MaterialPropertyName>("eT_beta"))),
   _h(getMaterialProperty<Real>(getParam<MaterialPropertyName>("h"))),
   _dh(getMaterialProperty<Real>("dh")),
   _a_val(declareProperty<Real>(getParam<MaterialPropertyName>("a"))),
   _da_de_val(declareProperty<Real>(getParam<MaterialPropertyName>("da_de"))),
   _da_dh_val(declareProperty<Real>(getParam<MaterialPropertyName>("da_dh")))
{
}

void
ScalarStrainJumpMaterial::computeQpProperties()
{//return the strain jump and its derivative wrt strain
  //a = (alpha -beta)
  _a_val[_qp] =    -(_mat_const_alpha[_qp]*(_grad_ux[_qp](0) - _eT_alpha[_qp]) - _mat_const_beta[_qp] * (_grad_ux[_qp](0) - _eT_beta[_qp]))/
                        (_mat_const_beta[_qp]* (1.0 - _h[_qp]) + _mat_const_alpha[_qp]*_h[_qp]);

  _da_de_val[_qp] = -(_mat_const_alpha[_qp] - _mat_const_beta[_qp])/
                        (_mat_const_beta[_qp]* (1.0 - _h[_qp]) + _mat_const_alpha[_qp]*_h[_qp]);

 _da_dh_val[_qp] =  _a_val[_qp]*(_mat_const_alpha[_qp] - _mat_const_beta[_qp])/
                        (_mat_const_beta[_qp]* (1.0 - _h[_qp]) + _mat_const_alpha[_qp]*_h[_qp]);
}
