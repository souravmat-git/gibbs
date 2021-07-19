//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
/* written by S.Chatterjee*/

#include "DerivativeStrainJump_3P_V1.h"
registerMooseObject("gibbsApp", DerivativeStrainJump_3P_V1);

template <>
InputParameters
validParams<DerivativeStrainJump_3P_V1>()
{
  InputParameters params = validParams<Material>();
  //params.addRequiredParam<MaterialPropertyName>("ds_de","Derivative wrt strain");
  //params.addRequiredParam<MaterialPropertyName>("da_dphi", "Derivative wrt phi");
  //params.addRequiredCoupledVar("eta", "phase-field variable");
  return params;
}

DerivativeStrainJump_3P_V1::DerivativeStrainJump_3P_V1(const InputParameters & parameters)
  : Material(parameters),
   //_grad_eta(coupledGradient("eta")),
   _n1(getMaterialProperty<RealVectorValue>("n1")),
   _n2(getMaterialProperty<RealVectorValue>("n2")),
   //Strain jumps
   _avec_alpha_beta(getMaterialProperty<RealVectorValue>("avec_alpha_beta")),
   _avec_beta_gamma(getMaterialProperty<RealVectorValue>("avec_beta_gamma")),
   //Obtain the prerequisite tensors
   //_lambda_star(getMaterialProperty<RankTwoTensor>("lambda_star")),
   //_L_hash(getMaterialProperty<RankTwoTensor>("L_hash")),
   _inv_D(getMaterialProperty<RankTwoTensor>("inv_D")),
   _S(getMaterialProperty<RankTwoTensor>("S")),
   //Stiffness of alpha, beta and gamma phases
   _alpha_stiffness(getMaterialProperty<RankFourTensor>("alpha_elasticity_tensor")),
   _beta_stiffness(getMaterialProperty<RankFourTensor>("beta_elasticity_tensor")),
   _gamma_stiffness(getMaterialProperty<RankFourTensor>("gamma_elasticity_tensor")),
   //Derivative of D with respect to phi_alpha, phi_beta and phi_gamma
   _dD_dphialpha(getMaterialProperty<RankTwoTensor>("dD_dphialpha")),
   _dD_dphibeta(getMaterialProperty<RankTwoTensor>("dD_dphibeta")),
   _dD_dphigamma(getMaterialProperty<RankTwoTensor>("dD_dphigamma")),
   //Derivative of strain jump with respect to strain to be calculated
   _ds_alpha_beta_de_val(declareProperty<RankFourTensor>("ds_alpha_beta_de")),
   _ds_beta_gamma_de_val(declareProperty<RankFourTensor>("ds_beta_gamma_de")),
   //Derivative of a_alpha_beta with respect to phi_alpha, phi_beta and phi_gamma
   _davec_alpha_beta_dphialpha_val(declareProperty<RealVectorValue>("davec_alpha_beta_dphialpha")),
   _davec_alpha_beta_dphibeta_val(declareProperty<RealVectorValue>("davec_alpha_beta_dphibeta")),
   _davec_alpha_beta_dphigamma_val(declareProperty<RealVectorValue>("davec_alpha_beta_dphigamma")),
   //Derivative of a_beta_gamma with respect to phi_alpha, phi_beta and phi_gamma
   _davec_beta_gamma_dphialpha_val(declareProperty<RealVectorValue>("davec_beta_gamma_dphialpha")),
   _davec_beta_gamma_dphibeta_val(declareProperty<RealVectorValue>("davec_beta_gamma_dphibeta")),
   _davec_beta_gamma_dphigamma_val(declareProperty<RealVectorValue>("davec_beta_gamma_dphigamma"))
{
}

void
DerivativeStrainJump_3P_V1::computeQpProperties()
{

  //const Real _norm_sq = _grad_eta[_qp].norm_sq();

  //if (_norm_sq < std::pow(libMesh::TOLERANCE, 3.0)){

  //  _ds_de_val[_qp].zero();
  //  _da_dphi_val[_qp].zero();

  //}
  //else{

    /*Declare a unit normal*/
  //   RealVectorValue _n  = -_grad_eta[_qp]/std::sqrt(_norm_sq);

  /*Declare two rank four tensors to store the difference*/
  RankFourTensor _stiff_alpha_beta, _stiff_beta_gamma;

  _stiff_alpha_beta = (_alpha_stiffness[_qp] - _beta_stiffness[_qp]);
  _stiff_beta_gamma = (_beta_stiffness[_qp] - _gamma_stiffness[_qp]);

  /*Declare two rank three tensors*/
  RankThreeTensor _da_alpha_beta_de, _da_beta_gamma_de;

  /*Set the rank three tensor to zero*/
  MathUtils::mooseSetToZero(_da_alpha_beta_de);
  MathUtils::mooseSetToZero(_da_beta_gamma_de);

  /*k,r,s are assumed to be free indices while i, p and a are repeated*/
    for(unsigned int k = 0; k < LIBMESH_DIM; k++)
      for(unsigned int r = 0; r < LIBMESH_DIM; r++)
        for(unsigned int s = 0; s < LIBMESH_DIM; s++)
          for(unsigned int p = 0; p < LIBMESH_DIM; p++)
            for(unsigned int a = 0; a < LIBMESH_DIM; a++)
              for(unsigned int i = 0; i < LIBMESH_DIM; i++)
              {
               _da_alpha_beta_de(k,r,s) += -_inv_D[_qp](p,k)*(_stiff_alpha_beta(p,a,r,s)*_n1[_qp](a));
                    //We have removed this terms since lambda_star will be zero
                    //+ _lambda_star[_qp](p,r) * _S[_qp](r,i) * _stiff_beta_gamma(i,a,r,s)*_n2[_qp](a));
              }

  /*Similarly, j,r and s are assumed to be free indices while i, a and k are repated*/
  for(unsigned int j = 0; j < LIBMESH_DIM; j++)
    for(unsigned int r = 0; r < LIBMESH_DIM; r++)
      for(unsigned int s = 0; s < LIBMESH_DIM; s++)
        for(unsigned int i = 0; i < LIBMESH_DIM; i++)
          for(unsigned int a = 0; a < LIBMESH_DIM; a++)
            for(unsigned int k = 0; k < LIBMESH_DIM; k++)
            {
              _da_beta_gamma_de(j,r,s) += -_S[_qp](j,i)*(_stiff_beta_gamma(i,a,r,s)*_n2[_qp](a));
              //We have removed this terms since L_hash will be zero
                            //+_L_hash[_qp](i,k) * _da_alpha_beta_de(k,r,s));
            }

  /*Then, determine the fourth rank-tensor*/
     for (unsigned int r = 0; r < LIBMESH_DIM; r++)
       for (unsigned int s = 0; s < LIBMESH_DIM; s++)
         for (unsigned int p = 0; p < LIBMESH_DIM; p++)
           for (unsigned int q = 0; q < LIBMESH_DIM; q++)
           {
              _ds_alpha_beta_de_val[_qp](r,s,p,q) =
                    (_da_alpha_beta_de(r,p,q)*_n1[_qp](s) + _n1[_qp](r)*_da_alpha_beta_de(s,p,q))/2.0;

              _ds_beta_gamma_de_val[_qp](r,s,p,q) =
                          (_da_beta_gamma_de(r,p,q)*_n2[_qp](s) + _n2[_qp](r)*_da_beta_gamma_de(s,p,q))/2.0;
           }

 //Calculate a_alpha_beta derivative with respect to phi_alpha
 //_davec_alpha_beta_dphialpha_val[_qp]  = -_inv_D[_qp] * _dD_dphialpha[_qp] * _avec_alpha_beta[_qp];
 //_davec_alpha_beta_dphibeta_val[_qp]   = -_inv_D[_qp] * _dD_dphibeta[_qp]  * _avec_alpha_beta[_qp];
 //_davec_alpha_beta_dphigamma_val[_qp]  = -_inv_D[_qp] * _dD_dphigamma[_qp] * _avec_alpha_beta[_qp];

 //_davec_beta_gamma_dphialpha_val[_qp]  = -_inv_D[_qp] * _dD_dphialpha[_qp] * _avec_beta_gamma[_qp];
 //_davec_beta_gamma_dphibeta_val[_qp]   = -_inv_D[_qp] * _dD_dphibeta[_qp]  * _avec_beta_gamma[_qp];
 //_davec_beta_gamma_dphigamma_val[_qp]  = -_inv_D[_qp] * _dD_dphigamma[_qp] * _avec_beta_gamma[_qp];

}
