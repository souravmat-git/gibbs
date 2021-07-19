//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "OnsagerFlux.h"

registerMooseObject("gibbsApp",OnsagerFlux);

template<>
InputParameters
validParams<OnsagerFlux>()
{
  InputParameters params = validParams<AuxKernel>();
  /**
    * Declare the options for a MooseEnum
    * Since the flux is a vector
    */
  MooseEnum component("x y z");
  //Use the MooseEnum to add a parameter called "component"
  params.addRequiredParam<MooseEnum>("component",component,"The desired component or direction.");
  params.addClassDescription("This is an aux kernel to calculate flux");
  //Add a "coupling field variable" to get a variable from input file
  params.addRequiredCoupledVar("B_diff_pot", "Diffusion Potential of component B");
  params.addRequiredCoupledVar("C_diff_pot", "Diffusion potential of component C");
  params.addRequiredParam<MaterialPropertyName>("diff_mobility_11", "Effect of mu1 on c1");
  params.addRequiredParam<MaterialPropertyName>("diff_mobility_12","Effect of mu2 on c1");
  
  return params;
}

OnsagerFlux::OnsagerFlux(const InputParameters & parameters)
  :AuxKernel(parameters),
  
  //This will automatically convert the MOOSEnum into an integer
  _component(getParam<MooseEnum>("component")),
  //Get the gradient of the variable
  _B_diff_pot_gradient(coupledGradient("B_diff_pot")),
  _C_diff_pot_gradient(coupledGradient("C_diff_pot")),
  _M_11(getMaterialProperty<Real>("diff_mobility_11")),
  _M_12(getMaterialProperty<Real>("diff_mobility_12"))
{
}
 
Real
OnsagerFlux::computeValue()
{
  //Access the gradient of the concentration
  // Then pull out the componnet of it we are looking for (x,y,z)
  //Note that getting a particular component of a gradient is done using the
  //paranthesis operator
  return -(_M_11[_qp]*(_B_diff_pot_gradient[_qp])(_component) + _M_12[_qp]*(_C_diff_pot_gradient[_qp])(_component));

}

