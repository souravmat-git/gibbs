//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//

#include "ComputeKtensor.h"
registerMooseObject("gibbsApp", ComputeKtensor);

//template <>
InputParameters
ComputeKtensor::validParams()
{
  InputParameters params = Material::validParams();
  //params.addRequiredParam<MaterialPropertyName>("stiffness_alpha",
  //                "Material constant for alpha phase");
  //params.addRequiredParam<MaterialPropertyName>("stiffness_beta",
  //                "Material constant for alpha phase");
  params.addRequiredParam<MaterialPropertyName>("K_name","K_tensor_name");
  return params;
}

ComputeKtensor::ComputeKtensor(const InputParameters & parameters)
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
ComputeKtensor::computeQpProperties(){

   _K_val[_qp].zero();
   _dK_dh_val[_qp].zero();

    for(unsigned int i=0; i< LIBMESH_DIM; ++i)
      for(unsigned int l=0; l<  LIBMESH_DIM; ++l)
        for(unsigned int k=0; k<  LIBMESH_DIM; ++k)
         for(unsigned int m=0; m<  LIBMESH_DIM; ++m){
            _K_val[_qp](i,l) += _n[_qp](k)*(_stiffness_alpha[_qp](k,i,l,m)*_h[_qp]
                      + _stiffness_beta[_qp](k,i,l,m)*(1.0 - _h[_qp])) *_n[_qp](m);

            _dK_dh_val[_qp](i,l) += _n[_qp](k)*(_stiffness_alpha[_qp](k,i,l,m)
                                       -_stiffness_beta[_qp](k,i,l,m))* _n[_qp](m);
         }

      //Check the determinant
      //std::cout << _K_val[_qp].det();
}
