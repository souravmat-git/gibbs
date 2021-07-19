//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SpecifyRadialDisplacementBC.h"
registerMooseObject("gibbsApp", SpecifyRadialDisplacementBC);

template <>
InputParameters
validParams<SpecifyRadialDisplacementBC>()
{
  InputParameters params = validParams<NodalBC>();
  params.addRequiredCoupledVar("disp_y", "Displacement component y");
  params.addRequiredParam<Real>("uR", "Radial displacement value");
  params.addRequiredParam<Point>("center_point",
                                 "Location of the center point of the cylindrical coordinates");
  return params;
}

SpecifyRadialDisplacementBC::SpecifyRadialDisplacementBC(const InputParameters & parameters)
  : NodalBC(parameters),
  _uy(coupledValue("disp_y")), // for the coupled variable
  _uy_var(coupled("disp_y")),
  _uR(getParam<Real>("uR")),
  _center_point(getParam<Point>("center_point"))
{
}

Real
SpecifyRadialDisplacementBC::computeQpResidual(){  

 Point _loc_from_center  = *_current_node - _center_point;
 //Returns the radial strain
 Real theta = std::atan2(_loc_from_center(1), _loc_from_center(0));

 //Based on coordinate transformation
 return _u[_qp]*std::cos(theta) + _uy[_qp]*std::sin(theta) - _uR;
}

Real
SpecifyRadialDisplacementBC::computeQpJacobian(){
 
  Point _loc_from_center  = *_current_node - _center_point;
  //Acts on the variable ux
  Real theta = std::atan2(_loc_from_center(1), _loc_from_center(0));
 
  return std::cos(theta);
}

Real
SpecifyRadialDisplacementBC::computeQpOffDiagJacobian(unsigned int jvar){

 Point _loc_from_center  = *_current_node - _center_point;
 //Acts on the variable uy
 Real theta = std::atan2(_loc_from_center(1), _loc_from_center(0));
  
 if(jvar == _uy_var){
    return std::sin(theta);
  }
  else
    return 0;  
}
