//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//

#include "ComputeInverseKtensor.h"
registerMooseObject("gibbsApp", ComputeInverseKtensor);

template <>
InputParameters
validParams<ComputeInverseKtensor>()
{
  InputParameters params = validParams<Material>();
  //params.addRequiredParam<MaterialPropertyName>("stiffness_alpha",
  //                "Material constant for alpha phase");
  //params.addRequiredParam<MaterialPropertyName>("stiffness_beta",
  //                "Material constant for alpha phase");
  params.addRequiredCoupledVar("eta","phase-field variable");
  params.addRequiredParam<MaterialPropertyName>("inv_K_name","inv_K_tensor_name");
  return params;
}

ComputeInverseKtensor::ComputeInverseKtensor(const InputParameters & parameters)
  : Material(parameters),
    _grad_eta(coupledGradient("eta")),
    //_n(getMaterialProperty<RealVectorValue>("n")),
   //Stiffness of alpha and beta phase
   _stiffness_alpha(getMaterialProperty<RankFourTensor>("alpha_elasticity_tensor")),
   _stiffness_beta(getMaterialProperty<RankFourTensor>("beta_elasticity_tensor")),
   //Interpolation function
   _h(getMaterialProperty<Real>("h")),
   //Compute the following properties
   _inv_K_val(declareProperty<RankTwoTensor>(getParam<MaterialPropertyName>("inv_K_name"))),
   _dK_dh_val(declareProperty<RankTwoTensor>("dK_dh"))
{
}

void
ComputeInverseKtensor::computeQpProperties(){

  const Real _norm_sq = _grad_eta[_qp].norm_sq();

  if (_norm_sq < std::pow(libMesh::TOLERANCE, 3.0) ){

	_inv_K_val[_qp].zero();
        _dK_dh_val[_qp].zero();

  }
  else {

   const RealVectorValue _n = - _grad_eta[_qp]/std::sqrt(_norm_sq);

   RankTwoTensor _K_val;

   _K_val.zero();
   _dK_dh_val[_qp].zero();

    for(unsigned int i=0; i< LIBMESH_DIM; ++i)
      for(unsigned int l=0; l<  LIBMESH_DIM; ++l)
        for(unsigned int k=0; k<  LIBMESH_DIM; ++k)
         for(unsigned int m=0; m<  LIBMESH_DIM; ++m){
            _K_val(i,l) += _n(k)*(_stiffness_alpha[_qp](k,i,l,m)*_h[_qp]
                      + _stiffness_beta[_qp](k,i,l,m)*(1.0 - _h[_qp])) *_n(m);

            _dK_dh_val[_qp](i,l) += _n(k)*(_stiffness_alpha[_qp](k,i,l,m)
                                       -_stiffness_beta[_qp](k,i,l,m))* _n(m);
         }


       _inv_K_val[_qp] = _K_val.inverse();

  }
}
