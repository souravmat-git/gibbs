//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DerivativePrerequisiteTensors_V1.h"
#include <fstream>
using namespace std;

registerMooseObject("gibbsApp", DerivativePrerequisiteTensors_V1);

template <>
InputParameters
validParams<DerivativePrerequisiteTensors_V1>()
{
  InputParameters params = validParams<Material>();
  params.addClassDescription("This material calculates the prerequsite"
                              "second-rank tensors for multiphase case");
  return params;
}

DerivativePrerequisiteTensors_V1::DerivativePrerequisiteTensors_V1(const InputParameters & parameters)
  : Material(parameters),
  _n1(getMaterialProperty<RealVectorValue>("n1")),
  _n2(getMaterialProperty<RealVectorValue>("n2")),
  //Stiffness of alpha and beta phase
  _alpha_elasticity_tensor(getMaterialProperty<RankFourTensor>("alpha_elasticity_tensor")),
  _beta_elasticity_tensor(getMaterialProperty<RankFourTensor>("beta_elasticity_tensor")),
  _gamma_elasticity_tensor(getMaterialProperty<RankFourTensor>("gamma_elasticity_tensor")),
  //Required derivatives of interpolation function h_alpha
  _dhalpha_dphialpha(getMaterialProperty<Real>("dhalpha_dphialpha")),
  _dhalpha_dphibeta(getMaterialProperty<Real>("dhalpha_dphibeta")),
  _dhalpha_dphigamma(getMaterialProperty<Real>("dhalpha_dphigamma")),
  //Required derivatives of interpolation function h_gamma
  _dhgamma_dphialpha(getMaterialProperty<Real>("dhgamma_dphialpha")),
  _dhgamma_dphibeta(getMaterialProperty<Real>("dhgamma_dphibeta")),
  _dhgamma_dphigamma(getMaterialProperty<Real>("dhgamma_dphigamma")),
  //Compute the following properties
  _lambda_hash(getMaterialProperty<RankTwoTensor>("lambda_hash")),
  //_lambda_star(getMaterialProperty<RankTwoTensor>("lambda_star")),
  //_L_hash(getMaterialProperty<RankTwoTensor>("L_hash")),
  _S(getMaterialProperty<RankTwoTensor>("S")),
  //Calculate
  _dD_dphialpha(declareProperty<RankTwoTensor>("dD_dphialpha")),
  _dD_dphibeta(declareProperty<RankTwoTensor>("dD_dphibeta")),
  _dD_dphigamma(declareProperty<RankTwoTensor>("dD_dphigamma"))
{
}

void
DerivativePrerequisiteTensors_V1::computeQpProperties(){

 //This class returns the prerequsite derivatives second rank tensors for 3P systems
 //for a given set of elastic tensors and unit normals
 RankFourTensor diff_alpha_beta, diff_beta_gamma;

 diff_alpha_beta = _alpha_elasticity_tensor[_qp] - _beta_elasticity_tensor[_qp];
 diff_beta_gamma = _beta_elasticity_tensor[_qp]  - _gamma_elasticity_tensor[_qp];

  //Derivative of lambda1 with respect to phi_alpha
 RankFourTensor dlambda1_dphialpha, dlambda1_dphibeta, dlambda1_dphigamma,
                dlambda2_dphialpha, dlambda2_dphibeta, dlambda2_dphigamma,
                dM1_dphialpha, dM1_dphibeta, dM1_dphigamma,
                dM2_dphialpha, dM2_dphibeta, dM2_dphigamma;

 //Derivative of lambda1 with respect to alpha, beta and gamma
 dlambda1_dphialpha  = -_dhalpha_dphialpha[_qp] * diff_alpha_beta;
 dlambda1_dphibeta   = -_dhalpha_dphibeta[_qp]  * diff_alpha_beta;
 dlambda1_dphigamma  = -_dhalpha_dphigamma[_qp] * diff_alpha_beta;

 //Derivative of lambda2 with respect to alpha, beta and gamma
 dlambda2_dphialpha  = _dhgamma_dphialpha[_qp] * diff_alpha_beta;
 dlambda2_dphibeta   = _dhgamma_dphibeta[_qp]  * diff_alpha_beta;
 dlambda2_dphigamma  = _dhgamma_dphigamma[_qp] * diff_alpha_beta;

 //Derivative of M1 with respect to alpha, beta and gamma
 dM1_dphialpha      = -_dhalpha_dphialpha[_qp] * diff_beta_gamma;
 dM1_dphibeta       = -_dhalpha_dphibeta[_qp]  * diff_beta_gamma;
 dM1_dphigamma      = -_dhalpha_dphigamma[_qp] * diff_beta_gamma;

 //Derivative of M2 with respect to alpha, beta and gamma
 dM2_dphialpha      = _dhgamma_dphialpha[_qp] * diff_beta_gamma;
 dM2_dphibeta       = _dhgamma_dphibeta[_qp]  * diff_beta_gamma;
 dM2_dphigamma      = _dhgamma_dphigamma[_qp] * diff_beta_gamma;

 RankTwoTensor _dlambda_hash_dphialpha, _dlambda_hash_dphibeta, _dlambda_hash_dphigamma,
               _dlambda_star_dphialpha, _dlambda_star_dphibeta, _dlambda_star_dphigamma,
               _dL_hash_dphialpha,  _dL_hash_dphibeta,  _dL_hash_dphigamma,
               _dL_star_dphialpha,  _dL_star_dphibeta,  _dL_star_dphigamma;

//Initialize the derivative of prerequsite tensors to zero
_dlambda_hash_dphialpha.zero();
_dlambda_hash_dphibeta.zero();
_dlambda_hash_dphigamma.zero();

/*
//lambda_star with respect to phi_alpha, phi_beta, phi_gamma
_dlambda_star_dphialpha.zero();
_dlambda_star_dphibeta.zero();
_dlambda_star_dphigamma.zero();
//L_hash with respect to phi_alpha, phi_beta, phi_gamma
_dL_hash_dphialpha.zero();
_dL_hash_dphibeta.zero();
_dL_hash_dphigamma.zero();

*/

//L_star with respect to phi_alpha, phi_beta and phi_gamma
_dL_star_dphialpha.zero();
_dL_star_dphibeta.zero();
_dL_star_dphigamma.zero();

 //Next, use for loops to calculate the tensors
 for(unsigned int i=0; i<LIBMESH_DIM; ++i)
   for(unsigned int k=0; k<LIBMESH_DIM; ++k)
      for(unsigned int l=0; l<LIBMESH_DIM; ++l)
        for(unsigned int j=0; j<LIBMESH_DIM; ++j)
        {

              _dlambda_hash_dphialpha(i,k) += _n1[_qp](l) * dlambda1_dphialpha(i,j,k,l) * _n1[_qp](j);
              _dlambda_hash_dphibeta(i,k)  += _n1[_qp](l) * dlambda1_dphibeta(i,j,k,l)  * _n1[_qp](j);
              _dlambda_hash_dphigamma(i,k) += _n1[_qp](l) * dlambda1_dphigamma(i,j,k,l) * _n1[_qp](j);
               /*
              _dlambda_star_dphialpha(i,k) += _n2[_qp](l) * dlambda2_dphialpha(i,j,k,l) * _n1[_qp](j);
              _dlambda_star_dphibeta(i,k)  += _n2[_qp](l) * dlambda2_dphibeta(i,j,k,l)  * _n1[_qp](j);
              _dlambda_star_dphigamma(i,k) += _n2[_qp](l) * dlambda2_dphigamma(i,j,k,l) * _n1[_qp](j);

              _dL_hash_dphialpha(i,k) += _n1[_qp](l) * dM1_dphialpha(i,j,k,l) * _n2[_qp](j);
              _dL_hash_dphibeta(i,k)  += _n1[_qp](l) * dM1_dphibeta(i,j,k,l)  * _n2[_qp](j);
              _dL_hash_dphigamma(i,k) += _n1[_qp](l) * dM1_dphigamma(i,j,k,l) * _n2[_qp](j);
              */
              _dL_star_dphialpha(i,k) += _n2[_qp](l) * dM2_dphialpha(i,j,k,l) * _n2[_qp](j);
              _dL_star_dphibeta(i,k)  += _n2[_qp](l) * dM2_dphibeta(i,j,k,l)  * _n2[_qp](j);
              _dL_star_dphigamma(i,k) += _n2[_qp](l) * dM2_dphigamma(i,j,k,l) * _n2[_qp](j);
        }

  RankTwoTensor _dS_dphialpha, _dS_dphibeta, _dS_dphigamma;
  //Calculate dS_dphialpha
  if (_dS_dphialpha.det()!=0){
  _dS_dphialpha  = _dL_star_dphialpha.inverse();
  }
  else{
   _dS_dphialpha.zero();
   }
  //dS_dphibeta
   if (_dS_dphibeta.det()!=0){
   _dS_dphibeta  = _dL_star_dphibeta.inverse();
   }
   else{
    _dS_dphibeta.zero();
    }
  //dS_dphigamma
    if (_dS_dphigamma.det()!=0){
    _dS_dphigamma  = _dL_star_dphigamma.inverse();
    }
    else{
     _dS_dphigamma.zero();
     }

  //std::cout << _dL_star_dphialpha.det();

  //print the value
  // if (_t_step)
  //{
  //  ofstream  file;
  // file.open("_dL_star_dphialpha.csv");

  //  for (unsigned int i = 0; i <3; ++i){
  //	  for (unsigned int j = 0; j <3; ++j){
  //          file << setprecision(6) << scientific << _dL_star_dphialpha(i,j) << ",";
  //     }
  //   file << "\n";
  //  }
  //}

  //Calculate dD_phialpha, dD_dphibeta and dD_dphigamma
  _dD_dphialpha[_qp] = _dlambda_hash_dphialpha;
                      //- _dlambda_star_dphialpha * _S[_qp] * _L_hash[_qp]
                      //- _lambda_star[_qp] * _dS_dphialpha * _L_hash[_qp]
                      //- _lambda_star[_qp] * _S[_qp] * _dL_hash_dphialpha;

  _dD_dphibeta[_qp] = _dlambda_hash_dphibeta;
                    //  - _dlambda_star_dphibeta * _S[_qp] * _L_hash[_qp]
                    //  - _lambda_star[_qp] * _dS_dphibeta * _L_hash[_qp]
                    //  - _lambda_star[_qp] * _S[_qp] * _dL_hash_dphibeta;

  _dD_dphigamma[_qp] = _dlambda_hash_dphigamma;
                    //  - _dlambda_star_dphigamma * _S[_qp] * _L_hash[_qp]
                    //  - _lambda_star[_qp] * _dS_dphigamma * _L_hash[_qp]
                    //  - _lambda_star[_qp] * _S[_qp] * _dL_hash_dphigamma;

}
