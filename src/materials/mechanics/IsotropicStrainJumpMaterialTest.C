//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//

#include "IsotropicStrainJumpMaterialTest.h"
registerMooseObject("gibbsApp", IsotropicStrainJumpMaterialTest);

//template <>
InputParameters
IsotropicStrainJumpMaterialTest::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<MaterialPropertyName>("total_strain",
                  "this tensor is calculated using strain-displacement relation");
  params.addRequiredParam<MaterialPropertyName>("inv_K_name",
                   "Inverse K tensor name");
  //params.addRequiredParam<MaterialPropertyName>("eigen_beta",
  //               "Transformation strain for beta phase");
  //params.addRequiredParam<MaterialPropertyName>("stiffness_alpha",
  //                "Material constant for alpha phase");
  //params.addRequiredParam<MaterialPropertyName>("stiffness_beta",
  //                "Material constant for alpha phase");
  params.addRequiredParam<MaterialPropertyName>("a_name","Jump in strain");
  return params;
}

IsotropicStrainJumpMaterialTest::IsotropicStrainJumpMaterialTest(const InputParameters & parameters)
  : Material(parameters),
    _n(getMaterialProperty<RealVectorValue>("n")),
    //Total strain
    _total_strain(getMaterialProperty<RankTwoTensor>(getParam<MaterialPropertyName>("total_strain"))),
    //inverse of rank two tensor
   _inv_K(getMaterialProperty<RankTwoTensor>(getParam<MaterialPropertyName>("inv_K_name"))),
    //Eigenstrains of alpha and beta phases
   _eigen_alpha(getMaterialProperty<RankTwoTensor>("alpha_eigen")),
   _eigen_beta(getMaterialProperty<RankTwoTensor>("beta_eigen")),
   //Stiffness of alpha and beta phase
   _alpha_lambda(getMaterialProperty<Real>("alpha_lambda")),
   _beta_lambda(getMaterialProperty<Real>("beta_lambda")),
   _alpha_mu(getMaterialProperty<Real>("alpha_mu")),
   _beta_mu(getMaterialProperty<Real>("beta_mu")),
   //Compute the following properties
   _a_val(declareProperty<RealVectorValue>(getParam<MaterialPropertyName>("a_name")))
{
}

void
IsotropicStrainJumpMaterialTest::computeQpProperties(){

    //Define the following quantities
    Real diff_lambda = (_alpha_lambda[_qp] - _beta_lambda[_qp]);
    Real diff_mu     = (_alpha_mu[_qp] - _beta_mu[_qp]);

    //Calculate X
    RealVectorValue _X =
       ( diff_lambda*_total_strain[_qp].trace()*_n[_qp] + 2.0*diff_mu*_total_strain[_qp]*_n[_qp]
       -_alpha_lambda[_qp]*_eigen_alpha[_qp].trace()*_n[_qp] - 2.0*_alpha_mu[_qp]*_eigen_alpha[_qp]*_n[_qp]
       +_beta_lambda[_qp]*_eigen_beta[_qp].trace()*_n[_qp]   + 2.0*_beta_mu[_qp] *_eigen_beta[_qp]*_n[_qp]);

    //Calculate a
    _a_val[_qp]  = - _inv_K[_qp] * _X;
}
