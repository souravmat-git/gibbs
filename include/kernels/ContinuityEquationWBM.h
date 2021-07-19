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

class ContinuityEquationWBM;

template <>
InputParameters validParams<ContinuityEquationWBM>();

/**
 * Kernel to implement the continuity
 * equation for mass conservation
 * based on WBM model. Note that continuity equation differs
 * between models such as KKS and GP based on the derivation of mu
 * This model acts on the mole fraction variable
 */
class ContinuityEquationWBM : public Kernel
{
public:
  ContinuityEquationWBM(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:

    Real interpolated_tf() const;
    RealGradient grad_mu() const;

    const VariableGradient & _grad_eta;
    unsigned int _eta_var;
    
    //const VariableValue & _diff_pot;
    //unsigned int _diff_pot_var;
  
    const MaterialProperty<Real> & _h;
    const MaterialProperty<Real> & _dh;
    const MaterialProperty<Real> & _d2h;
    
     //Keep in mind the free energy are functions of 
    //composition and not phase compositions
    //Note the difference in the naming of the two variables 
    //dfalpha_dc and df_dcalpha
    
    const MaterialProperty<Real> & _dfalpha_dc;
    const MaterialProperty<Real> & _d2falpha_dc2;
    const MaterialProperty<Real> & _dfbeta_dc;
    const MaterialProperty<Real> & _d2fbeta_dc2;
  
    /// Isotropic diffusion mobility indepnendent of concentration
    const MaterialProperty<Real> & _M;
};
