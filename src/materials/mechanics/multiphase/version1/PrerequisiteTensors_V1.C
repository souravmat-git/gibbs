//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//* Written by S.Chatterjee

#include "PrerequisiteTensors_V1.h"
registerMooseObject("gibbsApp", PrerequisiteTensors_V1);

template <>
InputParameters
validParams<PrerequisiteTensors_V1>()
{
  InputParameters params = validParams<VariableUnitNormalBase>();
  params.addClassDescription("This material calculates the prerequsite"
                              "second-rank tensors for multiphase case");
  return params;
}

PrerequisiteTensors_V1::PrerequisiteTensors_V1(const InputParameters & parameters)
  : VariableUnitNormalBase(parameters),
  //_n1(getMaterialProperty<RealVectorValue>("n1")),
  //_n2(getMaterialProperty<RealVectorValue>("n2")),
  //Stiffness of alpha and beta phase
  _alpha_elasticity_tensor(getMaterialProperty<RankFourTensor>("alpha_elasticity_tensor")),
  _beta_elasticity_tensor(getMaterialProperty<RankFourTensor>("beta_elasticity_tensor")),
  _gamma_elasticity_tensor(getMaterialProperty<RankFourTensor>("gamma_elasticity_tensor")),
  //Interpolation function
  _h_alpha(getMaterialProperty<Real>("h_alpha")),
  _h_beta(getMaterialProperty<Real>("h_beta")),
  _h_gamma(getMaterialProperty<Real>("h_gamma")),
  //Compute the following properties
  _lambda_hash(declareProperty<RankTwoTensor>("lambda_hash")),
  //_lambda_star(declareProperty<RankTwoTensor>("lambda_star")),
  //_L_hash(declareProperty<RankTwoTensor>("L_hash")),
  //_L_star(declareProperty<RankTwoTensor>("L_star")),
  _inv_D(declareProperty<RankTwoTensor>("inv_D")),
  _S(declareProperty<RankTwoTensor>("S"))
{
}

void
PrerequisiteTensors_V1::computeQpProperties(){

  RealVectorValue _n1, _n2;
  RankTwoTensor _L_star;

  //This class returns the prerequsite second rank tensors for 3P systems
  //for a given set of elastic tensors and unit normals
  if (unalpha_norm_sq() < std::pow(libMesh::TOLERANCE,3.0) ||
      unbeta_norm_sq()  < std::pow(libMesh::TOLERANCE,3.0)){
      //In the alpha-beta regions only the inverse of D is needed
      _inv_D[_qp].zero();

  } else {
      //First, define a unit normal for alpha/beta
      //Then, define a rank four tensor
      //Calculate the second rank tensor which will be non-zero only within the alpha-beta regions
      //Determine its inverse

      _n1 = - _grad_phibeta[_qp]/std::sqrt(unbeta_norm_sq());

      RankFourTensor lambda1 = _h_beta[_qp]  * _alpha_elasticity_tensor[_qp]
                             + _h_alpha[_qp] * _beta_elasticity_tensor[_qp]
                             + _h_gamma[_qp] * _alpha_elasticity_tensor[_qp];

      //Initialize lambda hash to zero
      _lambda_hash[_qp].zero();

      //Next, use for loops to calculate the tensors
      for(unsigned int i=0; i<LIBMESH_DIM; ++i)
         for(unsigned int k=0; k<LIBMESH_DIM; ++k)
           for(unsigned int l=0; l<LIBMESH_DIM; ++l)
             for(unsigned int j=0; j<LIBMESH_DIM; ++j){
                   _lambda_hash[_qp](i,k) += _n1(l) * lambda1(i,j,k,l)  * _n1(j);
           }

      //Declare a rank-two tensor Note that  lambda_star is zero
      //RankTwoTensor _D = _lambda_hash[_qp] - _lambda_star[_qp] * _S[_qp] * _L_hash[_qp];
      _inv_D[_qp] = _lambda_hash[_qp].inverse();

  }

  /***************************************************************************/

  if  (ungamma_norm_sq() < std::pow(libMesh::TOLERANCE,3.0))
  {   //Only the tensor S is needed for the beta/gamma region
    _S[_qp].zero();
  }
  else
  {

      _n2 = - _grad_phigamma[_qp]/std::sqrt(ungamma_norm_sq());

      RankFourTensor tensor_M2 =  _h_gamma[_qp]  * _beta_elasticity_tensor[_qp]
                                 + _h_alpha[_qp] * _gamma_elasticity_tensor[_qp]
                                 + _h_beta[_qp]  * _gamma_elasticity_tensor[_qp];

      //Initialize L_star to zero
      _L_star.zero();

      //Next, use for loops to calculate the tensors
      for(unsigned int i=0; i<LIBMESH_DIM; ++i)
        for(unsigned int k=0; k<LIBMESH_DIM; ++k)
          for(unsigned int l=0; l<LIBMESH_DIM; ++l)
            for(unsigned int j=0; j<LIBMESH_DIM; ++j){
              _L_star(i,k) += _n2(l) * tensor_M2(i,j,k,l) * _n2(j);
      }

     _S[_qp] = _L_star.inverse();
  }

  /* These two tensors are not needed
  RankFourTensor lambda2  =  _h_gamma[_qp] * _alpha_elasticity_tensor[_qp]
                           - _h_gamma[_qp] * _beta_elasticity_tensor[_qp];

  RankFourTensor tensor_M1 = _h_alpha[_qp] * _gamma_elasticity_tensor[_qp]
                           - _h_alpha[_qp] * _beta_elasticity_tensor[_qp];
  */

}
