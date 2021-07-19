//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PhaseFieldVelocity.h"

registerMooseObject("gibbsApp",PhaseFieldVelocity);

template<>
InputParameters
validParams<PhaseFieldVelocity>()
{
  InputParameters params = validParams<AuxKernel>();
  MooseEnum component("x y z");
  params.addRequiredParam<MooseEnum>("component",component,"The desired component or direction.");
  params.addRequiredCoupledVar("eta", "Phase field varible");
  params.addClassDescription("This is an aux kernel to calculate the phase field velcoity");  
  return params;
}

PhaseFieldVelocity::PhaseFieldVelocity(const InputParameters & parameters)
  :AuxKernel(parameters),
  _component(getParam<MooseEnum>("component")),
  _eta(coupledValue("eta")),
  _dot_eta(coupledDot("eta")),
  _grad_eta(coupledGradient("eta"))
{
}
 
Real
PhaseFieldVelocity::computeValue()
{
  //Access the gradient of the concentration
  // Then pull out the componnet of it we are looking for (x,y,z)
  //Note that getting a particular component of a gradient is done using the
  //paranthesis operator
  if (_eta[_qp] > 0.1 && _eta[_qp] < 0.9)
    return -(_dot_eta[_qp]/(_grad_eta[_qp](_component)));
  else
    return 0.0;

}

