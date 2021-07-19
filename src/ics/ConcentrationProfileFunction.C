//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.htmls

#include "ConcentrationProfileFunction.h"
#include "FEProblem.h"
#include "MooseMesh.h"

registerMooseObject("gibbsApp", ConcentrationProfileFunction);

template <>
InputParameters validParams<ConcentrationProfileFunction>()
{
  InputParameters params = validParams<EquilibriumProfileFunction>();
  params.addParam<Real>("xB_alpha", 0.1, "phase composition in alpha phase");
  params.addParam<Real>("xB_beta", 0.9, "phase composition in beta phase");
  return params;
}
ConcentrationProfileFunction::ConcentrationProfileFunction(const InputParameters & parameters)
  : EquilibriumProfileFunction(parameters),
   _xB_alpha(getParam<Real>("xB_alpha")),
   _xB_beta(getParam<Real>("xB_beta"))
{
}

Real
ConcentrationProfileFunction::h(const Point &p)
{
  const Real _factor = (std::sqrt(_W_height)/(std::sqrt(2.0)*_kappa));
  
  const Real _phi = 0.5*(1.0 - std::tanh(p(0)*_factor));
  
  //variable holding the value of interplation  
  const Real _h = std::pow(_phi,3.0)*(6.0*std::pow(_phi,2.0) - 15.0*_phi + 10);

  return (_h);
}

Real
ConcentrationProfileFunction::value(const Point &p)
{  
   //Note: unused parameters are commented, else it is a warning!
  //return : = c_beta*(1-h) + c_alpha*h
  return ( ConcentrationProfileFunction::h(p(0))*_xB_beta 
         +(1.0-ConcentrationProfileFunction::h(p(0)))*_xB_alpha); 
}
