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

class QMomentumBalance3DFun;
class Function;

template <>
InputParameters validParams<QMomentumBalance3DFun>();

/**
 * This kernel acts on the displacement variable
 * 1DKernel to implement momentum balance in two phase
 * material
 */
class QMomentumBalance3DFun : public Kernel
{
public:
  QMomentumBalance3DFun(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;


private:

  //The phase-field as a function
  const Function & _eta;

  //Global stress
  const MaterialProperty<RankTwoTensor> & _stress;
  
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
