//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//

#include "IsotropicDerivativeStrainJumpMaterial.h"
registerMooseObject("gibbsApp", IsotropicDerivativeStrainJumpMaterial);

template <>
InputParameters
validParams<IsotropicDerivativeStrainJumpMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<MaterialPropertyName>("ds_de","Derivative wrt strain");
  params.addRequiredParam<MaterialPropertyName>("da_dphi", "Derivative wrt phi");
  return params;
}

IsotropicDerivativeStrainJumpMaterial::IsotropicDerivativeStrainJumpMaterial(const InputParameters & parameters)
  : Material(parameters),
   _n(getMaterialProperty<RealVectorValue>("n")),
    //Jump in strain
    _a(getMaterialProperty<RealVectorValue>("a")),
   //Stiffness of alpha and beta phase
   _alpha_lambda(getMaterialProperty<Real>("alpha_lambda")),
   _beta_lambda(getMaterialProperty<Real>("beta_lambda")),
   _alpha_mu(getMaterialProperty<Real>("alpha_mu")),
   _beta_mu(getMaterialProperty<Real>("beta_mu")),
   //K tensor
   _inv_K(getMaterialProperty<RankTwoTensor>("inv_K")),
   //interpolation function
   _dh(getMaterialProperty<Real>("dh")),
   //Compute
   _ds_de_val(declareProperty<RankFourTensor>
                   (getParam<MaterialPropertyName>("ds_de"))),
   _da_dphi_val(declareProperty<RealVectorValue>(getParam<MaterialPropertyName>("da_dphi")))
{
}



void
IsotropicDerivativeStrainJumpMaterial::computeQpProperties()
{   
   //First, calculate a third rank tensor   
  const Real diff_lambda = (_alpha_lambda[_qp] - _beta_lambda[_qp]);
  const Real diff_mu     = (_alpha_mu[_qp] - _beta_mu[_qp]);
   
  RankTwoTensor I(RankTwoTensor::initIdentity);
  
  RankThreeTensor _da_de;
  MathUtils::mooseSetToZero(_da_de);
   
   //l,p,q are assumed to be free indices and i is repeated
  for (unsigned int l = 0; l < 3; l++)
    for (unsigned int p = 0; p < 3; p++)
      for (unsigned int q = 0; q <3; q++)
        for (unsigned int i = 0; i <3; i++)
        {
          _da_de(l,p,q) += - _inv_K[_qp](l,i)*(diff_lambda*I(p,q)*_n[_qp](i)
                                          + diff_mu*(I(i,p)*_n[_qp](q) + _n[_qp](p)*I(i,q))); 
        }
  
   //Then, determine the fourth rank-tensor
  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int k = 0; k < 3; k++)
      for (unsigned int p = 0; p < 3; p++)
        for (unsigned int q = 0; q < 3; q++)
        {
          _ds_de_val[_qp](i,k,p,q) = (_da_de(i,p,q)*_n[_qp](k) + _n[_qp](i)*_da_de(k,p,q))/2.0;
        }    
  
  
   /************  Calculate da/dh ***************************/
  
  RankTwoTensor P;
  P.vectorOuterProduct(_n[_qp], _n[_qp]);
 
  RankTwoTensor dK_dh = (diff_lambda + diff_mu)*P + diff_mu*I*_n[_qp].norm_sq(); 
   
  _da_dphi_val[_qp] = - _inv_K[_qp]*dK_dh*_a[_qp]* _dh[_qp];
}

