//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
// Forward Declarations
class ComputeOverall3PStress;

//MOOSe includes
#include "Material.h"

//template <>
//InputParameters validParams<ComputeOverall3PStress>();

class ComputeOverall3PStress : public Material
{
public:
  ComputeOverall3PStress(const InputParameters & parameters);

  static InputParameters validParams();

  //This class calculates the jump in strain or the strain difference
  //and its derivative wrt to overall strain

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:
    const MaterialProperty<RankTwoTensor> & _alpha_stress;
    const MaterialProperty<RankTwoTensor> & _beta_stress;
    const MaterialProperty<RankTwoTensor> & _gamma_stress;

    //Stiffness tensor of alpha and beta phase
    const MaterialProperty<RankFourTensor> & _alpha_stiffness;
    const MaterialProperty<RankFourTensor> & _beta_stiffness;
    const MaterialProperty<RankFourTensor> & _gamma_stiffness;

    //Strain jump and its derivatives
    const MaterialProperty<RankTwoTensor>  & _strain_jump_alpha_beta;
    const MaterialProperty<RankTwoTensor>  & _strain_jump_beta_gamma;

    //Derivative of strain jumps with respect to strain
    const MaterialProperty<RankFourTensor> & _ds_alpha_beta_de;
    const MaterialProperty<RankFourTensor> & _ds_beta_gamma_de;

    //Derivatives of strain jumps with respect to phi_alpha
    const MaterialProperty<RankTwoTensor> & _ds_alpha_beta_dphialpha;
    const MaterialProperty<RankTwoTensor> & _ds_beta_gamma_dphialpha;

   //Derivatives of strain jumps with respect to phi_beta
   const MaterialProperty<RankTwoTensor> & _ds_alpha_beta_dphibeta;
   const MaterialProperty<RankTwoTensor> & _ds_beta_gamma_dphibeta;

   //Derivatives of strain jumps with respect to phi_beta
   const MaterialProperty<RankTwoTensor> & _ds_alpha_beta_dphigamma;
   const MaterialProperty<RankTwoTensor> & _ds_beta_gamma_dphigamma;

    //interpolation functions
    const MaterialProperty<Real> & _h_alpha;
    const MaterialProperty<Real> & _h_beta;
    const MaterialProperty<Real> & _h_gamma;

    //Derivatives of interpolation function with respect to phi_alpha
    const MaterialProperty<Real> & _dhalpha_dphialpha;
    const MaterialProperty<Real> & _dhbeta_dphialpha;
    const MaterialProperty<Real> & _dhgamma_dphialpha;

    //Derivatives of interpolation function with respect to phi_beta
    const MaterialProperty<Real> & _dhalpha_dphibeta;
    const MaterialProperty<Real> & _dhgamma_dphibeta;

    //Derivatives of interpolation function with respect to gamma
    const MaterialProperty<Real> & _dhalpha_dphigamma;
    const MaterialProperty<Real> & _dhbeta_dphigamma;
    const MaterialProperty<Real> & _dhgamma_dphigamma;

    MaterialProperty<RankTwoTensor>  & _sigma_val;
    MaterialProperty<RankFourTensor> & _dsigma_de_val;

    MaterialProperty<RankTwoTensor>  & _dstress_dphialpha;
    MaterialProperty<RankTwoTensor>  & _dstress_dphibeta;
    MaterialProperty<RankTwoTensor>  & _dstress_dphigamma;
};
