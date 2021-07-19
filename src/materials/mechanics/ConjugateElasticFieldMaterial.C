//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//

#include "ConjugateElasticFieldMaterial.h"
registerMooseObject("gibbsApp", ConjugateElasticFieldMaterial);

template <>
InputParameters
validParams<ConjugateElasticFieldMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredCoupledVar("stress_x", "Stress in the x-direction");
  params.addRequiredParam<MaterialPropertyName>("mat_const", "Material constant for this phase");
  params.addRequiredParam<MaterialPropertyName>("eT", "Transformation strain");
  params.addRequiredParam<MaterialPropertyName>("ex", "The compatible strain in the x-direction");
  params.addRequiredParam<MaterialPropertyName>("compliance", "Compliance in the x-direction");
  params.addRequiredParam<MaterialPropertyName>("comp_energy","Elastic strain energy density of the phase");
  return params;
}

ConjugateElasticFieldMaterial::ConjugateElasticFieldMaterial(const InputParameters & parameters)
  : Material(parameters),
   _stress_x(coupledValue("stress_x")),
   _eT(getMaterialProperty<Real>(getParam<MaterialPropertyName>("eT"))),
   _mat_const(getMaterialProperty<Real>(getParam<MaterialPropertyName>("mat_const"))),
   _ex_val(declareProperty<Real>(getParam<MaterialPropertyName>("ex"))),
   _compliance_val(declareProperty<Real>(getParam<MaterialPropertyName>("compliance"))),
   _comp_energy_val(declareProperty<Real>(getParam<MaterialPropertyName>("comp_energy")))
{
} 

void
ConjugateElasticFieldMaterial::computeQpProperties()
{
    //x = 0, y = 1, z=2 (Strain in the material)
    _ex_val[_qp] = (_stress_x[_qp]/_mat_const[_qp]) + _eT[_qp];
    
    //Second derivative of the complementary strain energy
    _compliance_val[_qp] = (1.0/_mat_const[_qp]);
    
    //complementary strain energy
    _comp_energy_val[_qp] = -0.5*(_stress_x[_qp]*_stress_x[_qp])/_mat_const[_qp] 
                                - _eT[_qp]*_stress_x[_qp];
     
}
