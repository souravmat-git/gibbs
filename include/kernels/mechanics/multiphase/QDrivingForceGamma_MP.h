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

class QDrivingForceGamma_MP;

template<>
InputParameters validParams<QDrivingForceGamma_MP>();

/**
  *This class enforces the following
  * driving traction based on Aberyatne and Knowles
  **/

class QDrivingForceGamma_MP :  public Kernel
{
  public:
    QDrivingForceGamma_MP(const InputParameters & parameters);

  protected:

  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  private:

  Real driving_force() const;
  Real derivative_df() const;

  const VariableValue & _phase_beta;
  unsigned int _phase_beta_var;

  const VariableValue & _phase_alpha;
  unsigned int _phase_alpha_var;

  const MaterialProperty<Real> & _dhalpha_dphigamma;
  const MaterialProperty<Real> & _dhbeta_dphigamma;
  const MaterialProperty<Real> & _dhgamma_dphigamma;

  //Their second derivatives with resect to phigamma
  const MaterialProperty<Real> & _d2halpha_dphigamma2;
  const MaterialProperty<Real> & _d2hbeta_dphigamma2;

  //Second derivative with respect to phibeta
  const MaterialProperty<Real> & _d2halpha_dphigamma_dphibeta;
  const MaterialProperty<Real> & _d2hbeta_dphigamma_dphibeta;
  const MaterialProperty<Real> & _d2hgamma_dphigamma_dphibeta;

  //Second derivative with resect to phi_alpha
  const MaterialProperty<Real> & _d2halpha_dphigamma_dphialpha;
  const MaterialProperty<Real> & _d2hbeta_dphigamma_dphialpha;
  const MaterialProperty<Real> & _d2hgamma_dphigamma_dphialpha;

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
