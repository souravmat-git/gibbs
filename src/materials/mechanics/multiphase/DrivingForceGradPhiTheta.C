//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DrivingForceGradPhiTheta.h"
registerMooseObject("gibbsApp", DrivingForceGradPhiTheta);

template <>
InputParameters
validParams<DrivingForceGradPhiTheta>()
{
  InputParameters params = validParams<VariableUnitNormalBase>();
  return params;
}

DrivingForceGradPhiTheta::DrivingForceGradPhiTheta(const InputParameters & parameters)
  : VariableUnitNormalBase(parameters),
   _alpha_stress(getMaterialProperty<RankTwoTensor>("alpha_stress")),
   _beta_stress(getMaterialProperty<RankTwoTensor>("beta_stress")),
   _gamma_stress(getMaterialProperty<RankTwoTensor>("gamma_stress")),
   //Interpolation functions
   _h_alpha(getMaterialProperty<Real>("h_alpha")),
   _h_beta(getMaterialProperty<Real>("h_beta")),
   _h_gamma(getMaterialProperty<Real>("h_gamma")),
   //Jumpvectors
   _avec_alpha_beta(getMaterialProperty<RealVectorValue>("avec_alpha_beta")),
   _avec_beta_gamma(getMaterialProperty<RealVectorValue>("avec_beta_gamma")),
   //for alpha
   _alpha_dfbulk_dgradphi(declareProperty<RealVectorValue>("alpha_dfbulk_dgradphi")),
   _alpha_d2fbulk_dgradphi_dphi(declareProperty<RealVectorValue>("alpha_d2fbulk_dgradphi_dphi")),
   //for beta
   _beta_dfbulk_dgradphi(declareProperty<RealVectorValue>("beta_dfbulk_dgradphi")),
   _beta_d2fbulk_dgradphi_dphi(declareProperty<RealVectorValue>("beta_d2fbulk_dgradphi_dphi")),
   //for gamma
   _gamma_dfbulk_dgradphi(declareProperty<RealVectorValue>("gamma_dfbulk_dgradphi")),
   _gamma_d2fbulk_dgradphi_dphi(declareProperty<RealVectorValue>("gamma_d2fbulk_dgradphi_dphi"))
{
}

void
DrivingForceGradPhiTheta::computeQpProperties()
{
     //Initialize to zeros
    _alpha_dfbulk_dgradphi[_qp].zero();
    _alpha_d2fbulk_dgradphi_dphi[_qp].zero();

    //calculation for beta and gamma

    const RankTwoTensor I(RankTwoTensor::initIdentity);

    //This class returns the prerequsite second rank tensors for 3P systems
    //for a given set of elastic tensors and unit normals
    if (unalpha_norm_sq() < std::pow(libMesh::TOLERANCE,3.0) ||
        unbeta_norm_sq()  < std::pow(libMesh::TOLERANCE,3.0)){

        //In the bulk alpha-beta regions this contribution should be zero
        _beta_dfbulk_dgradphi[_qp].zero();
        _beta_d2fbulk_dgradphi_dphi[_qp].zero();

    }else {

      const Real _norm_cubebeta = std::pow(_grad_phibeta[_qp].norm(),3.0);

      //First define OG_beta
      RankTwoTensor OG_beta;
      OG_beta.vectorOuterProduct(_grad_phibeta[_qp], _grad_phibeta[_qp]);

      //Define the derivative of unit normal with respect to phase-field gradient
      RankTwoTensor _dn1_dgradphi_beta(-I/std::sqrt(unbeta_norm_sq())
                                       + OG_beta/_norm_cubebeta);

       // Then define
       //rem_ji = stress_jk * n_(k,i)
       //dL/dphi = M
       RankTwoTensor  _rem_beta = (_h_beta[_qp] * _h_alpha[_qp] * (_alpha_stress[_qp] - _beta_stress[_qp])
                                +  _h_gamma[_qp] * _h_alpha[_qp] * (_alpha_stress[_qp] - _gamma_stress[_qp]))
                                *  _dn1_dgradphi_beta;

       //RankTwoTensor  M = _dh[_qp]* (_alpha_stiffness[_qp] - _beta_stiffness[_qp])*_strain_jump[_qp] * _dn_dgradphi
        //                + (_h[_qp] * _alpha_stiffness[_qp] + (1.0 - _h[_qp])* _beta_stiffness[_qp])
        //                *_dstrainjump_dphi[_qp] *  _dn_dgradphi;


      //Initialize to zeros
      _beta_dfbulk_dgradphi[_qp].zero();
      _beta_d2fbulk_dgradphi_dphi[_qp].zero();

      for(unsigned int i=0; i<3; i++){
        for(unsigned int j=0; j<3; j++){
            _beta_dfbulk_dgradphi[_qp](i) +=  _rem_beta(j,i) * _avec_alpha_beta[_qp](j);
            //_beta_d2fbulk_dgradphi_dphi[_qp](i) +=  _dh[_qp] * (1.0 - 2.0 * _h[_qp])* L(j,i)*_a[_qp](j)
            //                              + _h[_qp] * (1.0 - _h[_qp]) * M(j,i)*_a[_qp](j)
            //                              + _h[_qp] * (1.0 - _h[_qp]) * L(j,i)* _da_dphi[_qp](j);
        }//end for
      }//end for
    }//end if

    //This class returns the prerequsite second rank tensors for 3P systems
    //for a given set of elastic tensors and unit normals
    if (ungamma_norm_sq() < std::pow(libMesh::TOLERANCE,3.0)){

        //In the bulk alpha-beta regions this contribution should be zero
        _gamma_dfbulk_dgradphi[_qp].zero();
        _gamma_d2fbulk_dgradphi_dphi[_qp].zero();

    } else {

      const Real _norm_cubegamma = std::pow(_grad_phigamma[_qp].norm(),3.0);

      //First define OG_beta
      RankTwoTensor OG_gamma;
      OG_gamma.vectorOuterProduct(_grad_phigamma[_qp], _grad_phigamma[_qp]);

      //Define the derivative of unit normal with respect to phase-field gradient
      RankTwoTensor _dn1_dgradphi_gamma( -I/std::sqrt(ungamma_norm_sq())
                                       + OG_gamma/_norm_cubegamma);

       // Then define
       //rem_ji = stress_jk * n_(k,i)
       //dL/dphi = M
       RankTwoTensor  _rem_gamma = (_h_gamma[_qp] * _h_alpha[_qp]* (_alpha_stress[_qp] - _gamma_stress[_qp])
                                  + _h_gamma[_qp]* _h_beta[_qp]  * (_beta_stress[_qp] - _gamma_stress[_qp]))
                                  * _dn1_dgradphi_gamma;

       //RankTwoTensor  M = _dh[_qp]* (_alpha_stiffness[_qp] - _beta_stiffness[_qp])*_strain_jump[_qp] * _dn_dgradphi
        //                + (_h[_qp] * _alpha_stiffness[_qp] + (1.0 - _h[_qp])* _beta_stiffness[_qp])
        //                *_dstrainjump_dphi[_qp] *  _dn_dgradphi;


      //Initialize to zeros
      _gamma_dfbulk_dgradphi[_qp].zero();
      _gamma_d2fbulk_dgradphi_dphi[_qp].zero();

      for(unsigned int i=0; i<3; i++){
        for(unsigned int j=0; j<3; j++){
            _gamma_dfbulk_dgradphi[_qp](i) +=  _rem_gamma(j,i) * _avec_beta_gamma[_qp](j);
            //_beta_d2fbulk_dgradphi_dphi[_qp](i) +=  _dh[_qp] * (1.0 - 2.0 * _h[_qp])* L(j,i)*_a[_qp](j)
            //                              + _h[_qp] * (1.0 - _h[_qp]) * M(j,i)*_a[_qp](j)
            //                              + _h[_qp] * (1.0 - _h[_qp]) * L(j,i)* _da_dphi[_qp](j);
        }//end for
      }//end for

  }//end if

}
