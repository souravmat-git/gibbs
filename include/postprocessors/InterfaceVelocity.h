/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#pragma once

#include "ChangeOverTimePostprocessor.h"
// Forward Declarations
class InterfaceVelocity;

template <>
InputParameters validParams<InterfaceVelocity>();

class InterfaceVelocity : public ChangeOverTimePostprocessor
{
public:
  InterfaceVelocity(const InputParameters & parameters);

protected:

  virtual Real getValue() override;
  Real _change_in_position;
  Real _v0;
};


