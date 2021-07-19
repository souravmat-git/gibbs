//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PrerequisiteVectors_V1.h"
registerMooseObject("gibbsApp", PrerequisiteVectors_V1);

template <>
InputParameters
validParams<PrerequisiteVectors_V1>()
{
  InputParameters params = validParams<VariableUnitNormalBase>();
  //params.addRequiredCoupledVar("displacements", "u v w");
  params.addClassDescription("This material calculates the prerequsite"
                              "second-rank tensors for multiphase case");
  return params;
}

PrerequisiteVectors_V1::PrerequisiteVectors_V1(const InputParameters & parameters)
  : VariableUnitNormalBase(parameters),
  //total strain
  _total_strain(getMaterialProperty<RankTwoTensor>("total_strain")),
  //unit normals
  _n1(getMaterialProperty<RealVectorValue>("n1")),
  _n2(getMaterialProperty<RealVectorValue>("n2")),
  //Stiffness of alpha, beta and gamma phases
  _alpha_elasticity_tensor(getMaterialProperty<RankFourTensor>("alpha_elasticity_tensor")),
  _beta_elasticity_tensor(getMaterialProperty<RankFourTensor>("beta_elasticity_tensor")),
  _gamma_elasticity_tensor(getMaterialProperty<RankFourTensor>("gamma_elasticity_tensor")),
  //Compute the following properties
  _m_alpha1(declareProperty<RealVectorValue>("m_alpha1")),
  _m_beta1(declareProperty<RealVectorValue>("m_beta1")),
  _psi_beta2(declareProperty<RealVectorValue>("psi_beta2")),
  _psi_gamma2(declareProperty<RealVectorValue>("psi_gamma2"))
{
}

void
PrerequisiteVectors_V1::computeQpProperties()
{
  /*
  if (_un1_norm_sq < std::pow(libMesh::TOLERANCE,3.0)){
     //||
      //_un2_norm_sq < std::pow(libMesh::TOLERANCE,3.0))
      //Set all required vectors to zero
      _m_alpha1[_qp].zero();
      _m_beta1[_qp].zero();
      _psi_beta2[_qp].zero();
      _psi_gamma2[_qp].zero();

  } else {

  */

    //n1 is a unit normal directed from beta to alpha phase
    //while n2 is a unit normal directed from gamma to beta phase
    //RealVectorValue _n1 = - _grad_phibeta[_qp]/std::sqrt(_un1_norm_sq);
    //RealVectorValue _n2 = - _grad_phigamma[_qp]/std::sqrt(_un2_norm_sq);

    _m_alpha1[_qp]   = _alpha_elasticity_tensor[_qp] * _total_strain[_qp] * _n1[_qp];
    _m_beta1[_qp]    = _beta_elasticity_tensor[_qp]  * _total_strain[_qp] * _n1[_qp];
    _psi_beta2[_qp]  = _beta_elasticity_tensor[_qp]  * _total_strain[_qp] * _n2[_qp];
    _psi_gamma2[_qp] = _gamma_elasticity_tensor[_qp] * _total_strain[_qp] * _n2[_qp];


  //}

}
