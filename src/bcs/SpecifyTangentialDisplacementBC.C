//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SpecifyTangentialDisplacementBC.h"
registerMooseObject("gibbsApp", SpecifyTangentialDisplacementBC);

template <>
InputParameters
validParams<SpecifyTangentialDisplacementBC>()
{
  InputParameters params = validParams<NodalBC>();
  params.addRequiredCoupledVar("disp_x", "Displacement component x");
  params.addRequiredParam<Real>("uT", "Specify tangential displacement");
  params.addRequiredParam<Point>("center_point",
                                 "Location of the center point of the cylindrical coordinates");
  return params;
}

SpecifyTangentialDisplacementBC::SpecifyTangentialDisplacementBC(const InputParameters & parameters)
  : NodalBC(parameters),
  _ux(coupledValue("disp_x")), // for the coupled variable
  _ux_var(coupled("disp_x")),
  _uT(getParam<Real>("uT")),
  _center_point(getParam<Point>("center_point"))
{
}

Real
SpecifyTangentialDisplacementBC::computeQpResidual(){

 Point _loc_from_center  = *_current_node - _center_point;

 Real theta = std::atan2(_loc_from_center(1), _loc_from_center(0));

  //Based on coordinate transformation
 return _ux[_qp]*std::cos(theta) + _u[_qp]*std::sin(theta) - _uT;
}

Real
SpecifyTangentialDisplacementBC::computeQpJacobian(){

  //Point _loc_from_center  = *_current_node - _center_point;
  //Returns the radial strain
  //Real theta = std::atan2(_loc_from_center(1), _loc_from_center(0));

  return 1;
}

Real
SpecifyTangentialDisplacementBC::computeQpOffDiagJacobian(unsigned int jvar){

 //Point _loc_from_center  = *_current_node - _center_point;
 //Returns the radial strain
 //Real theta = std::atan2(_loc_from_center(1), _loc_from_center(0));

 if(jvar == _ux_var){
    return 0;
  }
  else
    return 0;
}
