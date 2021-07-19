//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
#include "ACBulk.h"

class GMultiPhaseDoubleWell;

template <>
InputParameters validParams<GMultiPhaseDoubleWell>();

/**
 * Kernel to implement double well
 */
class GMultiPhaseDoubleWell : public ACBulk<Real>
{
public:
  GMultiPhaseDoubleWell(const InputParameters & parameters);

  protected:
    virtual Real computeDFDOP(PFFunctionType type) override;
    virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
    const VariableValue & _eta2;
    unsigned int _eta2_var;
    
    const VariableValue & _eta3;
    unsigned int _eta3_var;
  
    const VariableValue & _eta4;
    unsigned int _eta4_var;
  
    const VariableValue & _eta5;
    unsigned int _eta5_var;

    /// Mobility & Barrier height & gamma
    const MaterialProperty<Real> & _gamma;
    const MaterialProperty<Real> & _BH;

};
