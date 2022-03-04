//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DrivingForceGradPhi.h"
registerMooseObject("gibbsApp", DrivingForceGradPhi);

//template <>
InputParameters
DrivingForceGradPhi::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredCoupledVar("eta", "phase-field variable");
  params.addRequiredParam<MaterialPropertyName>("strain_jump_name",
                                            "Name of strain jump material");
  return params;
}

DrivingForceGradPhi::DrivingForceGradPhi(const InputParameters & parameters)
  : Material(parameters),
   _alpha_stress(getMaterialProperty<RankTwoTensor>("alpha_stress")),
   _beta_stress(getMaterialProperty<RankTwoTensor>("beta_stress")),
   //magnitude of a
   _a(getMaterialProperty<RealVectorValue>("a")),
   _grad_eta(coupledGradient("eta")),
   //_dn_dgradphi(getMaterialProperty<RankTwoTensor>("dn_dgradphi")),
   //Interpolation function
   _h(getMaterialProperty<Real>("h")),
   _dh(getMaterialProperty<Real>("dh")),
   _da_dphi(getMaterialProperty<RealVectorValue>("da_dphi")),
   _strain_jump(getMaterialProperty<RankTwoTensor>
                      (getParam<MaterialPropertyName>("strain_jump_name"))),
   //and its derivatives wrt to strain and phi
   _dstrainjump_dphi(getMaterialProperty<RankTwoTensor>("dstrainjump_dphi")),
   //Stiffness of alpha and beta phase
   _alpha_stiffness(getMaterialProperty<RankFourTensor>("alpha_elasticity_tensor")),
   _beta_stiffness(getMaterialProperty<RankFourTensor>("beta_elasticity_tensor")),
   //Compute the following properties
   _df_dgradphi(declareProperty<RealVectorValue>("df_dgradphi")),
   _d2f_dgradphi_dphi(declareProperty<RealVectorValue>("d2f_dgradphi_dphi"))
{
}

void
DrivingForceGradPhi::computeQpProperties()
{
   const Real _norm_sq = _grad_eta[_qp].norm_sq();
   const Real _norm_cube = std::pow(_grad_eta[_qp].norm(),3.0);
   const RankTwoTensor I(RankTwoTensor::initIdentity);

   if (_norm_sq < std::pow(libMesh::TOLERANCE,3.0)){

      _df_dgradphi[_qp].zero();
      _d2f_dgradphi_dphi[_qp].zero();

   } else {

    //First define OG
    RankTwoTensor OG;
    OG.vectorOuterProduct(_grad_eta[_qp], _grad_eta[_qp]);

    //Define another second rank tensor
    RankTwoTensor _dn_dgradphi;

    _dn_dgradphi = -I/std::sqrt(_norm_sq) + OG/_norm_cube;

    // Then define
    //L_ji = stress_jk * n_(k,i)
    //dL/dphi = M
    RankTwoTensor  L = (_alpha_stress[_qp] - _beta_stress[_qp])* _dn_dgradphi;
    RankTwoTensor  M = _dh[_qp]* (_alpha_stiffness[_qp] - _beta_stiffness[_qp])*_strain_jump[_qp] * _dn_dgradphi
                     + (_h[_qp] * _alpha_stiffness[_qp] + (1.0 - _h[_qp])* _beta_stiffness[_qp])
                                *_dstrainjump_dphi[_qp] *  _dn_dgradphi;

     //Initialize to zeros
    _df_dgradphi[_qp].zero();
    _d2f_dgradphi_dphi[_qp].zero();

    for(unsigned int i=0; i<3; i++){
      for(unsigned int j=0; j<3; j++){
        _df_dgradphi[_qp](i) += _h[_qp] * (1.0 - _h[_qp]) * L(j,i)*_a[_qp](j);
        _d2f_dgradphi_dphi[_qp](i) +=  _dh[_qp] * (1.0 - 2.0 * _h[_qp])* L(j,i)*_a[_qp](j)
                                     + _h[_qp] * (1.0 - _h[_qp]) * M(j,i)*_a[_qp](j)
                                     + _h[_qp] * (1.0 - _h[_qp]) * L(j,i)* _da_dphi[_qp](j);
      }
    }

  }
}
