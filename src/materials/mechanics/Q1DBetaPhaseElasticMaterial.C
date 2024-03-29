//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//

#include "Q1DBetaPhaseElasticMaterial.h"

registerMooseObject("gibbsApp", Q1DBetaPhaseElasticMaterial);

template <>
InputParameters
validParams<Q1DBetaPhaseElasticMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredCoupledVar("ux", "Displacement variable");
  params.addRequiredParam<MaterialPropertyName>("a", "Jump in strain");
  params.addRequiredParam<MaterialPropertyName>("h", "Interpolation Function");
  params.addRequiredParam<MaterialPropertyName>("eT", "Transformation strain");
  params.addRequiredParam<MaterialPropertyName>("mat_const", "Material property of phase");
  params.addRequiredParam<MaterialPropertyName>("ex", "Strain in the x-direction");
  params.addRequiredParam<MaterialPropertyName>("sx", "Stress in the x-direction");
  params.addRequiredParam<MaterialPropertyName>("elastic_energy","Phase strain energy density");
  return params;
}

Q1DBetaPhaseElasticMaterial::Q1DBetaPhaseElasticMaterial(const InputParameters & parameters)
  : Material(parameters),
   _grad_ux(coupledGradient("ux")),
   _a(getMaterialProperty<Real>(getParam<MaterialPropertyName>("a"))),
   _h(getMaterialProperty<Real>(getParam<MaterialPropertyName>("h"))),
   _eT(getMaterialProperty<Real>(getParam<MaterialPropertyName>("eT"))),
   _mat_const(getMaterialProperty<Real>(getParam<MaterialPropertyName>("mat_const"))),
   _ex_val(declareProperty<Real>(getParam<MaterialPropertyName>("ex"))),
   _sx_val(declareProperty<Real>(getParam<MaterialPropertyName>("sx"))),
   _elastic_energy_val(declareProperty<Real>(getParam<MaterialPropertyName>("elastic_energy")))
{
} 

void
Q1DBetaPhaseElasticMaterial::computeQpProperties(){  

  //return the stress and the elastic strain energy
  //Here, ex is the compatible strain
   _ex_val[_qp] = (_grad_ux[_qp](0) + _a[_qp] * (1.0- _h[_qp]));
   _sx_val[_qp] = _mat_const[_qp]* (_ex_val[_qp]  - _eT[_qp]);
   _elastic_energy_val[_qp] = 0.5*_mat_const[_qp]*(_ex_val[_qp] -_eT[_qp])*(_ex_val[_qp] -_eT[_qp]);  
}

