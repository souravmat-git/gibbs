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
class EqualPhaseTractionVectorXFun;
class Function;

template <>
InputParameters validParams<Kernel>();

/**
 * Enforce the equality of traction vector within each phase
 * in the x-direction
 */
class EqualPhaseTractionVectorXFun : public Kernel
{
  public:
    EqualPhaseTractionVectorXFun(const InputParameters & parameters);

  protected:
  
    virtual Real computeQpResidual() override;
    virtual Real computeQpJacobian() override;
    virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
   

    /// coupled variable for grad_ux_{beta}
    const  VariableValue & _ex_beta;
    unsigned int _ex_beta_var;
  
    ///Stress within each phase
    const MaterialProperty<Real> & _sx_alpha;
    const MaterialProperty<Real> & _sx_beta;
    
    //Elastic modulus of each phase
    const MaterialProperty<Real> & _mat_const_alpha;
    const MaterialProperty<Real> & _mat_const_beta;
};
