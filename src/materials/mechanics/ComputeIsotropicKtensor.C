//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//

#include "ComputeIsotropicKtensor.h"
#include "RankTwoTensor.h"
registerMooseObject("gibbsApp", ComputeIsotropicKtensor);

template <>
InputParameters
validParams<ComputeIsotropicKtensor>()
{
  InputParameters params = validParams<Material>();
  //params.addRequiredParam<MaterialPropertyName>("stiffness_alpha", 
  //                "Material constant for alpha phase");
  //params.addRequiredParam<MaterialPropertyName>("stiffness_beta", 
  //                "Material constant for alpha phase");
  params.addRequiredParam<MaterialPropertyName>("K_name","K_tensor_name");
  return params;
}

ComputeIsotropicKtensor::ComputeIsotropicKtensor(const InputParameters & parameters)
  : Material(parameters),
    _n(getMaterialProperty<RealVectorValue>("n")),
   //Stiffness of alpha and beta phase
   _stiffness_alpha(getMaterialProperty<RankFourTensor>("alpha_elasticity_tensor")),
   _stiffness_beta(getMaterialProperty<RankFourTensor>("beta_elasticity_tensor")),
   //Interpolation function  
   _h(getMaterialProperty<Real>("h")),
   //Compute the following properties
   _K_val(declareProperty<RankTwoTensor>(getParam<MaterialPropertyName>("K_name"))),
   _dK_dh_val(declareProperty<RankTwoTensor>("dK_dh"))
{
} 

void
ComputeIsotropicKtensor::computeQpProperties()
{

    _K_val[_qp].zero();
    _dK_dh_val[_qp].zero();
    
    const RankTwoTensor I(RankTwoTensor::initIdentity);
    RankTwoTensor P;
    
    P.vectorOuterProduct(_n[_qp], _n[_qp]);
    
    //lambda = C_1122 = C_12
    const Real lambda_alpha  = _stiffness_alpha[_qp](0,0,1,1);
    const Real lambda_beta   = _stiffness_beta[_qp](0,0,1,1);
    
    //mu = C_2323 = C_44 (Voigt Notation)                 
    const Real mu_alpha = _stiffness_alpha[_qp](1,2,1,2);
    const Real mu_beta  = _stiffness_beta[_qp](1,2,1,2);  

   
    const Real lambda = _h[_qp]* lambda_alpha +  (1.0- _h[_qp]) * lambda_beta ;                     
    const Real mu     = _h[_qp]* mu_alpha +  (1.0 - _h[_qp]) * mu_beta;
                              
    const Real diff_lambda = (lambda_alpha - lambda_beta);
    const Real diff_mu     = (mu_alpha - mu_beta);
                              
                                                 
    _K_val[_qp] = (lambda + mu)* P + mu*I*_n[_qp].norm_sq();  
    _dK_dh_val[_qp] = (diff_lambda + diff_mu)* P + diff_mu*I*_n[_qp].norm_sq();
        
    //Check                    
    //Diagonal components 
    //_K_val[_qp](0,0) = (lambda + mu)*_n[_qp](0)*_n[_qp](0) + mu*_n[_qp].norm_sq();
    //_K_val[_qp](1,1) = (lambda + mu)*_n[_qp](1)*_n[_qp](1) + mu*_n[_qp].norm_sq();
    //_K_val[_qp](2,2) = (lambda + mu)*_n[_qp](2)*_n[_qp](2) + mu*_n[_qp].norm_sq();
    
    //off-diagonal components 
    //_K_val[_qp](0,1) = (lambda + mu)*_n[_qp](0)*_n[_qp](1);
    //_K_val[_qp](0,2) = (lambda + mu)*_n[_qp](0)*_n[_qp](2);
    //_K_val[_qp](1,2) = (lambda + mu)*_n[_qp](1)*_n[_qp](2);
    
    //Due to symmetry
    //_K_val[_qp](1,0) = _K_val[_qp](0,1);
    //_K_val[_qp](2,0) = _K_val[_qp](0,2);
    //_K_val[_qp](2,1) = _K_val[_qp](1,2);
    

}

