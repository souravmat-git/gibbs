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

class QMomentumBalance1D_3P;
template <>
InputParameters validParams<QMomentumBalance1D_3P>();

/**
 * This kernel acts on the displacement variable
 * 1DKernel to implement momentum balance in two phase
 * material
 */
class QMomentumBalance1D_3P : public Kernel
{
public:
  QMomentumBalance1D_3P(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;


private:

  //Real interp_modulus() const;

  //The phase-field variables
  //const VariableValue & _;
  //unsigned int _eta_var;

  //Stress in the beta phase and the stress in the alpha phase
  const MaterialProperty<Real> & _sx_alpha;

  //material constant (lambda + 2mu) in the beta phase and the alpha phase
  const MaterialProperty<Real> & _mat_const_alpha;

  //Strain jump and its derivative
  //const MaterialProperty<Real> & _a;
  //const MaterialProperty<Real> & _da_de;
  //const MaterialProperty<Real> & _da_dh;

  //Interpolation function
  //const MaterialProperty<Real> & _h;
  //const MaterialProperty<Real> & _dh;
};
