//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DerivativePrerequisiteTensors.h"
registerMooseObject("gibbsApp", DerivativePrerequisiteTensors);

//template <>
InputParameters
DerivativePrerequisiteTensors::validParams()
{
  InputParameters params = VariableUnitNormalBase::validParams();
  params.addClassDescription("This material calculates the prerequsite"
                              "second-rank tensors for multiphase case");
  return params;
}

DerivativePrerequisiteTensors::DerivativePrerequisiteTensors(const InputParameters & parameters)
  : VariableUnitNormalBase(parameters),
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
   //Calculate in the alpha-beta regions
   _dD_dphialpha(declareProperty<RankTwoTensor>("dD_dphialpha")),
   _dD_dphibeta(declareProperty<RankTwoTensor>("dD_dphibeta")),
   _dD_dphigamma(declareProperty<RankTwoTensor>("dD_dphigamma")),
   //Calculate in the beta-gamma regions
   _dS_dphialpha(declareProperty<RankTwoTensor>("dS_dphialpha")),
   _dS_dphibeta(declareProperty<RankTwoTensor>("dS_dphibeta")),
   _dS_dphigamma(declareProperty<RankTwoTensor>("dS_dphigamma"))
{
}

void
DerivativePrerequisiteTensors::computeQpProperties(){

  if (unalpha_norm_sq() < std::pow(libMesh::TOLERANCE,3.0) ||
      unbeta_norm_sq()  < std::pow(libMesh::TOLERANCE,3.0)){

        _dD_dphialpha[_qp].zero();
        _dD_dphibeta[_qp].zero();
        _dD_dphigamma[_qp].zero();

  } else {

    RealVectorValue _n1 = - _grad_phibeta[_qp]/std::sqrt(unbeta_norm_sq());

    RankFourTensor diff_alpha_beta =
           (_alpha_elasticity_tensor[_qp] - _beta_elasticity_tensor[_qp]);

    RankFourTensor dlambda1_dphialpha, dlambda1_dphibeta, dlambda1_dphigamma;

    //Derivative of lambda1 with respect to alpha, beta and gamma
    dlambda1_dphialpha  = -_dhalpha_dphialpha[_qp] * diff_alpha_beta;
    dlambda1_dphibeta   = -_dhalpha_dphibeta[_qp]  * diff_alpha_beta;
    dlambda1_dphigamma  = -_dhalpha_dphigamma[_qp] * diff_alpha_beta;

    //Initialize to zero
    _dD_dphialpha[_qp].zero();
    _dD_dphibeta[_qp].zero();
    _dD_dphigamma[_qp].zero();

    //Next, use for loops to calculate the tensors
    for(unsigned int i=0; i<LIBMESH_DIM; ++i)
      for(unsigned int k=0; k<LIBMESH_DIM; ++k)
         for(unsigned int l=0; l<LIBMESH_DIM; ++l)
           for(unsigned int j=0; j<LIBMESH_DIM; ++j){
                 _dD_dphialpha[_qp](i,k) += _n1(l) * dlambda1_dphialpha(i,j,k,l) * _n1(j);
                 _dD_dphibeta[_qp](i,k)  += _n1(l) * dlambda1_dphibeta(i,j,k,l)  * _n1(j);
                 _dD_dphigamma[_qp](i,k) += _n1(l) * dlambda1_dphigamma(i,j,k,l) * _n1(j);
           }//for loop ends
   }//if loop ends

  if (ungamma_norm_sq() < std::pow(libMesh::TOLERANCE,3.0)){
     _dS_dphialpha[_qp].zero();
     _dS_dphibeta[_qp].zero();
     _dS_dphigamma[_qp].zero();

  } else {

    RealVectorValue _n2 = - _grad_phigamma[_qp]/std::sqrt(ungamma_norm_sq());

    RankFourTensor diff_beta_gamma
                = _beta_elasticity_tensor[_qp]  - _gamma_elasticity_tensor[_qp];

    RankFourTensor dM2_dphialpha, dM2_dphibeta, dM2_dphigamma;

    //Derivative of M2 with respect to alpha, beta and gamma
    dM2_dphialpha      = _dhgamma_dphialpha[_qp] * diff_beta_gamma;
    dM2_dphibeta       = _dhgamma_dphibeta[_qp]  * diff_beta_gamma;
    dM2_dphigamma      = _dhgamma_dphigamma[_qp] * diff_beta_gamma;

    RankTwoTensor _dL_star_dphialpha, _dL_star_dphibeta, _dL_star_dphigamma;

    //Initialize _dL_star_dphitheta tensors to zero
    _dL_star_dphialpha.zero();
    _dL_star_dphibeta.zero();
    _dL_star_dphigamma.zero();

    //Next, use for loops to calculate the tensors
    for(unsigned int i=0; i<LIBMESH_DIM; ++i)
      for(unsigned int k=0; k<LIBMESH_DIM; ++k)
         for(unsigned int l=0; l<LIBMESH_DIM; ++l)
           for(unsigned int j=0; j<LIBMESH_DIM; ++j){
                 _dL_star_dphialpha(i,k) += _n2(l) * dM2_dphialpha(i,j,k,l) * _n2(j);
                 _dL_star_dphibeta(i,k)  += _n2(l) * dM2_dphibeta(i,j,k,l)  * _n2(j);
                 _dL_star_dphigamma(i,k) += _n2(l) * dM2_dphigamma(i,j,k,l) * _n2(j);
           }//end for loop

     //Note that dL_star can be zero due to the derivative
     // of interpolation function hence we need to check
     //Take the inverse by checking for the zero
     if (_dL_star_dphialpha.det()!= 0)
       _dS_dphialpha[_qp] = _dL_star_dphialpha.inverse();
     else
       _dS_dphialpha[_qp].zero();

     if (_dL_star_dphibeta.det()!= 0)
       _dS_dphibeta[_qp]  = _dL_star_dphibeta.inverse();
     else
       _dS_dphibeta[_qp].zero();

     if (_dL_star_dphigamma.det()!= 0)
      _dS_dphigamma[_qp] = _dL_star_dphigamma.inverse();
     else
      _dS_dphigamma[_qp].zero();

   }//end if
}
