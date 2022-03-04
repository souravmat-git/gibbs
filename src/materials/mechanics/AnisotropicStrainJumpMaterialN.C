//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AnisotropicStrainJumpMaterialN.h"
registerMooseObject("gibbsApp", AnisotropicStrainJumpMaterialN);

//template <>
InputParameters
AnisotropicStrainJumpMaterialN::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredCoupledVar("displacements", "u v w");
  //params.addRequiredParam<MaterialPropertyName>("eigen_alpha",
  //                "Transformation strain for alpha phase");
  //params.addRequiredParam<MaterialPropertyName>("eigen_beta",
  //               "Transformation strain for beta phase");
  //params.addRequiredParam<MaterialPropertyName>("stiffness_alpha",
  //                "Material constant for alpha phase");
  //params.addRequiredParam<MaterialPropertyName>("stiffness_beta",
  //                "Material constant for alpha phase");
  params.addRequiredCoupledVar("eta", "phase-field variable");
  params.addRequiredParam<MaterialPropertyName>("a_name","Jump in strain");
  return params;
}

AnisotropicStrainJumpMaterialN::AnisotropicStrainJumpMaterialN(const InputParameters & parameters)
  : Material(parameters),
   _grad_eta(coupledGradient("eta")),
   //_n(getMaterialProperty<RealVectorValue>("n")),
   _ndisp(coupledComponents("displacements")),
   //_disp(coupledValues("displacements")),
   //Vector array of displacement gradients
   _grad_disp(coupledGradients("displacements")),
   //Eigenstrains of alpha and beta phase
   _eigen_alpha(getMaterialProperty<RankTwoTensor>("alpha_eigen")),
   _eigen_beta(getMaterialProperty<RankTwoTensor>("beta_eigen")),
   //Stiffness of alpha and beta phase
   _stiffness_alpha(getMaterialProperty<RankFourTensor>("alpha_elasticity_tensor")),
   _stiffness_beta(getMaterialProperty<RankFourTensor>("beta_elasticity_tensor")),
   //inverse K tensor
   _inv_K(getMaterialProperty<RankTwoTensor>("inv_K")),
   //Compute the following properties
   _a_val(declareProperty<RealVectorValue>(getParam<MaterialPropertyName>("a_name")))
{
  //unsigned int ndisp = coupledComponents("displacements");
  //Check if number of displacements = mesh dimension
  if (_ndisp != _mesh.dimension())
      mooseError(
        "The number of variables supplied in 'disp' must match the mesh dimension.");

  //_disp.resize(3, &_zero);
  _grad_disp.resize(3, &_grad_zero);

  //Depending on the mesh dimension the components of disp gradient to zero
  //for(unsigned i = ndisp; i<3; i++)
  //  _grad_disp[i] = &_grad_zero;
}

void
AnisotropicStrainJumpMaterialN::computeQpProperties(){

   const Real _norm_sq = _grad_eta[_qp].norm_sq();

   if (_norm_sq < std::pow(libMesh::TOLERANCE, 3.0)){
        _a_val[_qp].zero();
   }
   else{

   //Define the unit normal
   RealVectorValue _n = - _grad_eta[_qp]/std::sqrt(_norm_sq);

   //Define the displacement gradient tensor which is in general not symmetric
   const auto _disp_grad_tensor = RankTwoTensor::initializeFromRows((*_grad_disp[0])[_qp],
                                    (*_grad_disp[1])[_qp],
                                    (*_grad_disp[2])[_qp]);

   //Assuming infinitesimal deformation, the total strain is given by
   RankTwoTensor _total_strain = (_disp_grad_tensor + _disp_grad_tensor.transpose()) / 2.0;


   RankTwoTensor M = (_stiffness_alpha[_qp] - _stiffness_beta[_qp]) * _total_strain
                    -(_stiffness_alpha[_qp] * _eigen_alpha[_qp] - _stiffness_beta[_qp] * _eigen_beta[_qp]);

   RealVectorValue _X = M.transpose() * _n;

   //X_i = [C_iklm]*\epsilon_lm*n_k
   //      - (C^{\alpha}_iklm*\eigen^{\alpha}_lm - C^{\beta}_{iklm}*\eigen^{\beta_lm})*n_k
   //RealVectorValue _X =
   //                  (_stiffness_alpha[_qp] - _stiffness_beta[_qp]) * _total_strain * _n[_qp]
   //                 -(_stiffness_alpha[_qp] * _eigen_alpha[_qp] - _stiffness_beta[_qp] * _eigen_beta[_qp])*_n[_qp];

    _a_val[_qp]  = -_inv_K[_qp] * _X;
   }
}
