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

class QMomentumBalance3D_MP;
template <>
InputParameters validParams<QMomentumBalance3D_MP>();

class QMomentumBalance3D_MP : public Kernel
{
public:
  QMomentumBalance3D_MP(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;


private:

  //The phase-field variables
  const VariableValue & _phase_alpha;
  const VariableValue & _phase_beta;
  const VariableValue & _phase_gamma;

  unsigned int _phase_alpha_var;
  unsigned int _phase_beta_var;
  unsigned int _phase_gamma_var;

  //Global stress
  const MaterialProperty<RankTwoTensor> & _stress;

  //Derivative with respect to phase alpha
  const MaterialProperty<RankTwoTensor> & _dstress_dphialpha;
  const MaterialProperty<RankTwoTensor> & _dstress_dphibeta;
  const MaterialProperty<RankTwoTensor> & _dstress_dphigamma;

  //Jacobian mult = d\sigma_ij/d\epsilon_ij
  const MaterialProperty<RankFourTensor> & _Jacobian_mult;

  //The direction of displacement variable that this kernel acts on
  //_component = 0 = x, _component = 1 = y, _component = 2 = z
  const unsigned int _component;

  //Number of coupled displacements  variable
  unsigned int _ndisp;

  //Displacement variables IDs
  std::vector<unsigned int> _disp_var;

};
