//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/* written by S.Chatterjee*/

#include "DerivativeStrainJump_3P.h"
registerMooseObject("gibbsApp", DerivativeStrainJump_3P);

//template <>
InputParameters
DerivativeStrainJump_3P::validParams()
{
  InputParameters params = VariableUnitNormalBase::validParams();
  return params;
}

DerivativeStrainJump_3P::DerivativeStrainJump_3P(const InputParameters & parameters)
  : VariableUnitNormalBase(parameters),
   //Strain jumps
   _avec_alpha_beta(getMaterialProperty<RealVectorValue>("avec_alpha_beta")),
   _psi_beta2(getMaterialProperty<RealVectorValue>("psi_beta2")),
   _psi_gamma2(getMaterialProperty<RealVectorValue>("psi_gamma2")),
   //Obtain the prerequisite tensors
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
   //Derivative of S with respect to phi_alpha, phi_beta and phi_gamma
   _dS_dphialpha(getMaterialProperty<RankTwoTensor>("dS_dphialpha")),
   _dS_dphibeta(getMaterialProperty<RankTwoTensor>("dS_dphibeta")),
   _dS_dphigamma(getMaterialProperty<RankTwoTensor>("dS_dphigamma")),
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
DerivativeStrainJump_3P::computeQpProperties()
{

  if (unalpha_norm_sq() < std::pow(libMesh::TOLERANCE,3.0) ||
      unbeta_norm_sq()  < std::pow(libMesh::TOLERANCE,3.0)){

     _ds_alpha_beta_de_val[_qp].zero();
     _davec_alpha_beta_dphialpha_val[_qp].zero();
     _davec_alpha_beta_dphibeta_val[_qp].zero();
     _davec_alpha_beta_dphigamma_val[_qp].zero();

  } else {

    RealVectorValue _n1 = -_grad_phibeta[_qp]/std::sqrt(unbeta_norm_sq());

    RankFourTensor _diff_stiff_alpha_beta =
                        (_alpha_stiffness[_qp] - _beta_stiffness[_qp]);

    /*Declare a rank three tensors*/
    RankThreeTensor _da_alpha_beta_de;

    /*Set the rank three tensor to zero*/
    MathUtils::mooseSetToZero(_da_alpha_beta_de);

    /*k,r,s are assumed to be free indices while i, p and a are repeated*/
      for(unsigned int k = 0; k < LIBMESH_DIM; k++)
        for(unsigned int r = 0; r < LIBMESH_DIM; r++)
          for(unsigned int s = 0; s < LIBMESH_DIM; s++)
            for(unsigned int p = 0; p < LIBMESH_DIM; p++)
              for(unsigned int a = 0; a < LIBMESH_DIM; a++){
                  _da_alpha_beta_de(k,r,s) += -_inv_D[_qp](p,k)*(_diff_stiff_alpha_beta(p,a,r,s)*_n1(a));
                }//end for

     //Then, determine the fourth rank-tensor*/
      for(unsigned int r = 0; r < LIBMESH_DIM; r++)
        for(unsigned int s = 0; s < LIBMESH_DIM; s++)
          for(unsigned int p = 0; p < LIBMESH_DIM; p++)
            for (unsigned int q = 0; q < LIBMESH_DIM; q++){
                    _ds_alpha_beta_de_val[_qp](r,s,p,q) =
                                  (_da_alpha_beta_de(r,p,q)*_n1(s)
                                +  _n1(r)*_da_alpha_beta_de(s,p,q))/2.0;
              }//end for

  //Finally, calculate a_alpha_beta derivative with respect to phi_alpha, phi_beta and phi_gamma
  _davec_alpha_beta_dphialpha_val[_qp]  = -_inv_D[_qp] * _dD_dphialpha[_qp] * _avec_alpha_beta[_qp];
  _davec_alpha_beta_dphibeta_val[_qp]   = -_inv_D[_qp] * _dD_dphibeta[_qp]  * _avec_alpha_beta[_qp];
  _davec_alpha_beta_dphigamma_val[_qp]  = -_inv_D[_qp] * _dD_dphigamma[_qp] * _avec_alpha_beta[_qp];

}//end if

if (ungamma_norm_sq() < std::pow(libMesh::TOLERANCE,3.0)){

    _ds_beta_gamma_de_val[_qp].zero();
    _davec_beta_gamma_dphialpha_val[_qp].zero();
    _davec_beta_gamma_dphibeta_val[_qp].zero();
    _davec_beta_gamma_dphigamma_val[_qp].zero();

} else {

  RealVectorValue _n2 = - _grad_phigamma[_qp]/std::sqrt(ungamma_norm_sq());

  /*Declare a rank four tensors to store the difference*/
  RankFourTensor _diff_stiff_beta_gamma =
                   (_beta_stiffness[_qp] - _gamma_stiffness[_qp]);

  /*Declare a rank-three tensor*/
  RankThreeTensor _da_beta_gamma_de;

  //Initialize the rank-three tensor to zero
  MathUtils::mooseSetToZero(_da_beta_gamma_de);

  for(unsigned int j = 0; j < LIBMESH_DIM; j++)
    for(unsigned int r = 0; r < LIBMESH_DIM; r++)
      for(unsigned int s = 0; s < LIBMESH_DIM; s++)
        for(unsigned int i = 0; i < LIBMESH_DIM; i++)
          for(unsigned int a = 0; a < LIBMESH_DIM; a++)
          {
              _da_beta_gamma_de(j,r,s) += -_S[_qp](j,i)*(_diff_stiff_beta_gamma(i,a,r,s)*_n2(a));
          }//end for

    /*Then, determine the fourth rank-tensor*/
  for(unsigned int r = 0; r < LIBMESH_DIM; r++)
    for(unsigned int s = 0; s < LIBMESH_DIM; s++)
      for(unsigned int p = 0; p < LIBMESH_DIM; p++)
        for(unsigned int q = 0; q < LIBMESH_DIM; q++)
        {
            _ds_beta_gamma_de_val[_qp](r,s,p,q) =
                        (_da_beta_gamma_de(r,p,q)*_n2(s)
                        + _n2(r)*_da_beta_gamma_de(s,p,q))/2.0;
        }//end for

   RealVectorValue delta_psi = (_psi_beta2[_qp] - _psi_gamma2[_qp]);

   //Calculate avec_beta_gamma derivative with resect to phi_alpha, phi_beta and phi_gamma
   _davec_beta_gamma_dphialpha_val[_qp]  = -_dS_dphialpha[_qp] * delta_psi;
   _davec_beta_gamma_dphibeta_val[_qp]   = -_dS_dphibeta[_qp]  * delta_psi;
   _davec_beta_gamma_dphigamma_val[_qp]  = -_dS_dphigamma[_qp] * delta_psi;

  }//end if
}//end the member function
