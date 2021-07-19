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
class TwoDimEqualPhaseTractionVectorYFun;

template <>
InputParameters validParams<Kernel>();

/**
 * Enforce the equality of traction vector within each phase
 * in the x-direction
 */
class TwoDimEqualPhaseTractionVectorYFun : public Kernel
{
  public:
    TwoDimEqualPhaseTractionVectorYFun(const InputParameters & parameters);

  protected:
  
    virtual Real computeQpResidual() override;
    virtual Real computeQpJacobian() override;
    virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

    /// coupled variable eyy_beta
    //const  VariableValue & _eyy_beta;
    unsigned int _eyy_beta_var;
    
    ///Stress xy of each phase
    const MaterialProperty<Real> & _syy_alpha;
    const MaterialProperty<Real> & _syy_beta;
    
    //C66 each phase
    const MaterialProperty<Real> & _C22_alpha;
    const MaterialProperty<Real> & _C22_beta;

};
