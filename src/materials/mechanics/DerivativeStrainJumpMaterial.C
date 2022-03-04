//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
/* written by S.Chatterjee*/

#include "DerivativeStrainJumpMaterial.h"
registerMooseObject("gibbsApp", DerivativeStrainJumpMaterial);

//template <>
InputParameters
DerivativeStrainJumpMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<MaterialPropertyName>("ds_de","Derivative wrt strain");
  params.addRequiredParam<MaterialPropertyName>("da_dphi", "Derivative wrt phi");
  params.addRequiredCoupledVar("eta", "phase-field variable");
  return params;
}

DerivativeStrainJumpMaterial::DerivativeStrainJumpMaterial(const InputParameters & parameters)
  : Material(parameters),
   _grad_eta(coupledGradient("eta")),
   /*_n(getMaterialProperty<RealVectorValue>("n")),*/
   //Jump in strain
    _a(getMaterialProperty<RealVectorValue>("a")),
   //Eigenstrains of alpha and beta phase
   _eigen_alpha(getMaterialProperty<RankTwoTensor>("alpha_eigen")),
   _eigen_beta(getMaterialProperty<RankTwoTensor>("beta_eigen")),
   //Stiffness of alpha and beta phase
   _stiffness_alpha(getMaterialProperty<RankFourTensor>("alpha_elasticity_tensor")),
   _stiffness_beta(getMaterialProperty<RankFourTensor>("beta_elasticity_tensor")),
   //K tensor
   _K(getMaterialProperty<RankTwoTensor>("K")),
   _dK_dh(getMaterialProperty<RankTwoTensor>("dK_dh")),
   /*Derivative of interpolation function*/
   _dh(getMaterialProperty<Real>("dh")),
   _ds_de_val(declareProperty<RankFourTensor>
                   (getParam<MaterialPropertyName>("ds_de"))),
   _da_dphi_val(declareProperty<RealVectorValue>(getParam<MaterialPropertyName>("da_dphi")))
{
}

void
DerivativeStrainJumpMaterial::computeQpProperties()
{

  const Real _norm_sq = _grad_eta[_qp].norm_sq();

  if (_norm_sq < std::pow(libMesh::TOLERANCE,2.0)){

    _ds_de_val[_qp].zero();
    _da_dphi_val[_qp].zero();

  }
  else{

    /*Declare a unit normal*/
     RealVectorValue _n  = -_grad_eta[_qp]/std::sqrt(_norm_sq);

    /*Declare a rank three tensor*/
    RankThreeTensor _da_de;

    /*Set the rank three tensor to zero*/
    MathUtils::mooseSetToZero(_da_de);

    /*Compute the inverse of K tensor*/
    const RankTwoTensor _inv_K = _K[_qp].inverse();

    /*l,p,q are assumed to be free indices while i and k are repeated*/
    for (unsigned int l = 0; l < LIBMESH_DIM; l++)
      for (unsigned int p = 0; p < LIBMESH_DIM; p++)
        for (unsigned int q = 0; q <LIBMESH_DIM; q++)
          for (unsigned int i = 0; i <LIBMESH_DIM; i++)
            for (unsigned int k = 0; k < LIBMESH_DIM; k++){
              _da_de(l,p,q) += - _inv_K(l,i)*(_stiffness_alpha[_qp](p,q,k,i)
                                          - _stiffness_beta[_qp](p,q,k,i))*_n(k);
            }

     /*Then, determine the fourth rank-tensor*/
     for (unsigned int r = 0; r < LIBMESH_DIM; r++)
       for (unsigned int s = 0; s < LIBMESH_DIM; s++)
         for (unsigned int p = 0; p < LIBMESH_DIM; p++)
           for (unsigned int q = 0; q < LIBMESH_DIM; q++){
              _ds_de_val[_qp](r,s,p,q) = (_da_de(r,p,q)*_n(s) + _n(r)*_da_de(s,p,q))/2.0;
           }

    /************  Calculate da/dphi ***************************/
     _da_dphi_val[_qp] = -_inv_K*_dK_dh[_qp]*_a[_qp]*_dh[_qp];
  }
}
