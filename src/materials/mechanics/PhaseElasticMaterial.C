//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//

#include "PhaseElasticMaterial.h"
registerMooseObject("gibbsApp", PhaseElasticMaterial);

template <>
InputParameters
validParams<PhaseElasticMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredCoupledVar("ex", "Compatible strain for an arbitrary phase");
  params.addRequiredParam<MaterialPropertyName>("eT", "Transformation strain");
  params.addRequiredParam<MaterialPropertyName>("mat_const", "Material property of phase");
  params.addRequiredParam<MaterialPropertyName>("sx", "Stress in the x-direction");
  params.addRequiredParam<MaterialPropertyName>("elastic_energy","Elastic strain energy density of the phase");
  return params;
}

PhaseElasticMaterial::PhaseElasticMaterial(const InputParameters & parameters)
  : Material(parameters),
   _ex(coupledValue("ex")),
   _eT(getMaterialProperty<Real>(getParam<MaterialPropertyName>("eT"))),
   _mat_const(getMaterialProperty<Real>(getParam<MaterialPropertyName>("mat_const"))),
   _sx_val(declareProperty<Real>(getParam<MaterialPropertyName>("sx"))),
   _elastic_energy_val(declareProperty<Real>(getParam<MaterialPropertyName>("elastic_energy")))
{
} 

void
PhaseElasticMaterial::computeQpProperties(){  

  //return the stress and the elastic strain energy
  //Here, ex is the elastic strain i,e, ex = ec-eT
   _sx_val[_qp] = _mat_const[_qp]* _ex[_qp];
   _elastic_energy_val[_qp] = 0.5*_mat_const[_qp]*_ex[_qp]*_ex[_qp];  
}
