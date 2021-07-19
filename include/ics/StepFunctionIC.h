//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef STEPFUNCTION_H
#define STEPFUNCTION_H

#include "InitialCondition.h"
#include <cmath>
using namespace std;

// Forward Declarations
class StepFunctionIC;

template <>
InputParameters validParams<StepFunctionIC>();

/**
 * Makes initial condition which creates a linear ramp of the given variable
 * on the x-axis with specified side values
 */
class StepFunctionIC : public InitialCondition
{
public:

  StepFunctionIC(const InputParameters & parameters);
  
  virtual Real value(const Point & p) override;

protected:
  /**
   * The value of the variable at a point.
   */
  const Real _xlength;
  const Real _ylength;
  const Real _phase_comp_alpha;
  const Real _phase_comp_beta;
  const Real _phase_comp_gamma;
  const Real _phase_len;
};

#endif // RAMPIC_H
