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

// Forward Declarations
class TwoDimEqualPhaseTractionVectorXFun;

template <>
InputParameters validParams<Kernel>();
/**
 * Enforce the equality of traction vector within each phase
 * in the x-direction
 */
class TwoDimEqualPhaseTractionVectorXFun : public Kernel
{
  public:
    TwoDimEqualPhaseTractionVectorXFun(const InputParameters & parameters);

  protected:
  
    virtual Real computeQpResidual() override;
    virtual Real computeQpJacobian() override;
    virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

    /// coupled variable exy_beta
    //const  VariableValue & _exy_beta;
    unsigned int _exy_beta_var;
    
    ///Stress xy of each phase
    const MaterialProperty<Real> & _sxy_alpha;
    const MaterialProperty<Real> & _sxy_beta;
    
    //C66 each phase
    const MaterialProperty<Real> & _C66_alpha;
    const MaterialProperty<Real> & _C66_beta;

};
