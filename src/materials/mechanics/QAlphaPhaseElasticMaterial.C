//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//

#include "QAlphaPhaseElasticMaterial.h"
registerMooseObject("gibbsApp", QAlphaPhaseElasticMaterial);

//template <>
InputParameters
QAlphaPhaseElasticMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<std::string>("base_name","The phase name");
  //params.addRequiredParam<MaterialPropertyName>("total_strain",
  //                "this tensor is calculated using strain-displacement relation");
  //params.addRequiredParam<MaterialPropertyName>("eigen_strain",
  //               "this tensor is generally a user-defined input eigenstrain name");
  //params.addRequiredParam<MaterialPropertyName>("elasticity_tensor",
  //                                               "Stiffness tensor");
  //params.addRequiredParam<MaterialPropertyName>("strain_jump", "Jump in strain");
  //params.addRequiredParam<MaterialPropertyName>("h", "Interpolation Function");
  return params;
}

QAlphaPhaseElasticMaterial::QAlphaPhaseElasticMaterial(const InputParameters & parameters)
  : Material(parameters),
   //Get the phase name which will be appended to the material properties
   _phase_name(getParam<std::string>("base_name")),
   //Total strain
   _total_strain(getMaterialProperty<RankTwoTensor>("total_strain")),
   //Eigenstrain,  strain jump of beta phase
   _alpha_eigen(getMaterialProperty<RankTwoTensor>("alpha_eigen")),
   _strain_jump(getMaterialProperty<RankTwoTensor>("strain_jump")),
   //elasticity tensor to calculate stress
   _alpha_stiffness(getMaterialProperty<RankFourTensor>("alpha_elasticity_tensor")),
    //interpolation function
   _h(getMaterialProperty<Real>("h")),
   //Following properties will be calculated by this material class
   _elastic_strain_val(declareProperty<RankTwoTensor>(_phase_name + "_elastic_strain")),
   _stress_val(declareProperty<RankTwoTensor>(_phase_name + "_stress")),
   //_Jacobian_mult(declareProperty<RankFourTensor>(_phase_name + "_Jacobian_mult")),
   _elastic_energy_val(declareProperty<Real>(_phase_name + "_strain_energy_density"))
{
}

void
QAlphaPhaseElasticMaterial::computeQpProperties(){

  _elastic_strain_val[_qp].zero();
  _stress_val[_qp].zero();

  //Alpha elastic strain = (total) compatible strain - alpha_eigenstrain
  //+  strain jump * interpolation
   _elastic_strain_val[_qp] = _total_strain[_qp]
                    + _h[_qp] * _strain_jump[_qp] - _alpha_eigen[_qp];

  //Following this calculate the phase stress
   _stress_val[_qp] = _alpha_stiffness[_qp]* _elastic_strain_val[_qp];

   //J_ijpq = C_ijpq + h*C_ijkl * ds_kl/de_pq, where ds_kl/de_pq = M_klpq
  //_Jacobian_mult[_qp] = _alpha_stiffness[_qp]
   //                     + _h[_qp]* _alpha_stiffness[_qp] * _ds_de[_qp];

  //Finally, calculate the elastic strain energy
   _elastic_energy_val[_qp] =
        _stress_val[_qp].doubleContraction((_elastic_strain_val)[_qp]) / 2.0;
}
