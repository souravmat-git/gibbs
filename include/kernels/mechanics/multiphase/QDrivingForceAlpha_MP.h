//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
#include "Kernel.h"

class QDrivingForceAlpha_MP;

template<>
InputParameters validParams<QDrivingForceAlpha_MP>();

/**
  *This class enforces the following
  * driving traction based on Aberyatne and Knowles
  **/

class QDrivingForceAlpha_MP :  public Kernel
{
  public:
    QDrivingForceAlpha_MP(const InputParameters & parameters);

  protected:

  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  private:

  Real driving_force() const;
  Real derivative_df() const;

  const VariableValue & _phase_beta;
  unsigned int _phase_beta_var;

  const VariableValue & _phase_gamma;
  unsigned int _phase_gamma_var;

  //Their first derivatives with resect to alpha
  const MaterialProperty<Real> & _dhalpha_dphialpha;
  const MaterialProperty<Real> & _dhbeta_dphialpha;
  const MaterialProperty<Real> & _dhgamma_dphialpha;

  //Their second derivatives with resect to phi_alpha
  const MaterialProperty<Real> & _d2hbeta_dphialpha2;
  const MaterialProperty<Real> & _d2hgamma_dphialpha2;

  //wrt beta
  const MaterialProperty<Real> & _d2halpha_dphialpha_dphibeta;
  const MaterialProperty<Real> & _d2hbeta_dphialpha_dphibeta;
  const MaterialProperty<Real> & _d2hgamma_dphialpha_dphibeta;

  //wrt gamma
  const MaterialProperty<Real> & _d2halpha_dphialpha_dphigamma;
  const MaterialProperty<Real> & _d2hbeta_dphialpha_dphigamma;
  const MaterialProperty<Real> & _d2hgamma_dphialpha_dphigamma;

  //Strain energies
  const MaterialProperty<Real> & _alpha_strain_energy;
  const MaterialProperty<Real> & _beta_strain_energy;
  const MaterialProperty<Real> & _gamma_strain_energy;

  //stress tensor
  const MaterialProperty<RankTwoTensor> & _stress;

  //Strain jumps
  const MaterialProperty<RankTwoTensor> & _strain_jump_alpha_beta;
  const MaterialProperty<RankTwoTensor> & _strain_jump_beta_gamma;

  //phase-field mobility assumed to be equal in all phases
  const MaterialProperty<Real> & _L;

  //nd_factor
  const MaterialProperty<Real> & _nd_factor;

};
