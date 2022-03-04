//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//

#include "AnisotropicStrainJumpMaterialTest.h"
registerMooseObject("gibbsApp", AnisotropicStrainJumpMaterialTest);

//template <>
InputParameters
AnisotropicStrainJumpMaterialTest::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<MaterialPropertyName>("total_strain",
                  "this tensor is calculated using strain-displacement relation");
  //params.addRequiredParam<MaterialPropertyName>("eigen_alpha",
  //                "Transformation strain for alpha phase");
  //params.addRequiredParam<MaterialPropertyName>("eigen_beta",
  //               "Transformation strain for beta phase");
  //params.addRequiredParam<MaterialPropertyName>("stiffness_alpha",
  //                "Material constant for alpha phase");
  //params.addRequiredParam<MaterialPropertyName>("stiffness_beta",
  //                "Material constant for alpha phase");
  params.addRequiredParam<MaterialPropertyName>("a_name","Jump in strain");
  return params;
}

AnisotropicStrainJumpMaterialTest::AnisotropicStrainJumpMaterialTest(const InputParameters & parameters)
  : Material(parameters),
    _n(getMaterialProperty<RealVectorValue>("n")),
    //Total strain
    _total_strain(getMaterialProperty<RankTwoTensor>(getParam<MaterialPropertyName>("total_strain"))),
   //Eigenstrains of alpha and beta phase
   _eigen_alpha(getMaterialProperty<RankTwoTensor>("alpha_eigen")),
   _eigen_beta(getMaterialProperty<RankTwoTensor>("beta_eigen")),
   //Stiffness of alpha and beta phase
   _stiffness_alpha(getMaterialProperty<RankFourTensor>("alpha_elasticity_tensor")),
   _stiffness_beta(getMaterialProperty<RankFourTensor>("beta_elasticity_tensor")),
   //K tensor
   _K(getMaterialProperty<RankTwoTensor>("K")),
   //Compute the following properties
   _a_val(declareProperty<RealVectorValue>(getParam<MaterialPropertyName>("a_name")))
{
}

void
AnisotropicStrainJumpMaterialTest::computeQpProperties()
{
  //X_i = [C_iklm]*\epsilon_lm*n_k
  //      - (C^{\alpha}_iklm*\eigen^{\alpha}_lm - C^{\beta}_{iklm}*\eigen^{\beta_lm})*n_k
  const RealVectorValue _X =
                  (_stiffness_alpha[_qp] - _stiffness_beta[_qp]) * _total_strain[_qp] * _n[_qp]
                 -(_stiffness_alpha[_qp] * _eigen_alpha[_qp] - _stiffness_beta[_qp] * _eigen_beta[_qp])*_n[_qp];

   _a_val[_qp]  = -_K[_qp].inverse() * _X;

}
