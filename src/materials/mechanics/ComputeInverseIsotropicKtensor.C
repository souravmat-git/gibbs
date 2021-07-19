//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//

#include "ComputeInverseIsotropicKtensor.h"
#include "RankTwoTensor.h"
registerMooseObject("gibbsApp", ComputeInverseIsotropicKtensor);

template <>
InputParameters
validParams<ComputeInverseIsotropicKtensor>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<MaterialPropertyName>("inv_K_name","inv_K_tensor_name");
  return params;
}

ComputeInverseIsotropicKtensor::ComputeInverseIsotropicKtensor(const InputParameters & parameters)
  : Material(parameters),
    _n(getMaterialProperty<RealVectorValue>("n")),
   //Stiffness of alpha and beta phase
   _alpha_lambda(getMaterialProperty<Real>("alpha_lambda")),
   _beta_lambda(getMaterialProperty<Real>("beta_lambda")),
   _alpha_mu(getMaterialProperty<Real>("alpha_mu")),
   _beta_mu(getMaterialProperty<Real>("beta_mu")),
   //Interpolation function  
   _h(getMaterialProperty<Real>("h")),
   //Compute the following properties
   _inv_K_val(declareProperty<RankTwoTensor>(getParam<MaterialPropertyName>("inv_K_name"))),
   _dK_dh_val(declareProperty<RankTwoTensor>("dK_dh"))
{
} 

void
ComputeInverseIsotropicKtensor::computeQpProperties(){

    //Define the following quantities
    Real lambda = _h[_qp]* _alpha_lambda[_qp] +  (1.0- _h[_qp]) * _beta_lambda[_qp] ;                     
    Real mu     = _h[_qp]* _alpha_mu[_qp] +  (1.0 - _h[_qp]) * _beta_mu[_qp];
                              
    Real det_K = mu*(lambda + 2.0*mu)*std::pow(_n[_qp].norm_sq(),2.0);
    
    RankTwoTensor I(RankTwoTensor::initIdentity);
    RankTwoTensor P;
    
    P.vectorOuterProduct(_n[_qp], _n[_qp]);
    
    _inv_K_val[_qp] = ((lambda + 2.0*mu)*I*_n[_qp].norm_sq() -(lambda + mu)*P)/det_K;
    
    _dK_dh_val[_qp].zero();
  
}

