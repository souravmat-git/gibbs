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
class TractionEquality1D;


template <>
InputParameters validParams<Kernel>();

/**
 * Enforce the equality of stress within each phase
 * in the x-direction
 */
class TractionEquality1D : public Kernel
{
  public:
    TractionEquality1D(const InputParameters & parameters);

  protected:
   
    virtual Real computeQpResidual() override;
    virtual Real computeQpJacobian() override;
    virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
    
    Real rev_interp_mat_const() const ;
    
    unsigned int _eta_var;
    unsigned int _ux_var;
   
    ///Stress within each phase
    const MaterialProperty<Real> & _sx_alpha;
    const MaterialProperty<Real> & _sx_beta;
    
    //Elastic modulus of each phase
    const MaterialProperty<Real> & _mat_const_alpha;
    const MaterialProperty<Real> & _mat_const_beta;
    
    //Interpolation material is required
    const MaterialProperty<Real> & _h;
    const MaterialProperty<Real> & _dh;
};
