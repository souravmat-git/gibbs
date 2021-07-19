//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CylindricalRankOneAux.h"
registerMooseObject("gibbsApp", CylindricalRankOneAux);

InputParameters
CylindricalRankOneAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Takes a nodal vector variable and outputs its component in cylindrical coordinates");
  params.addRequiredCoupledVar("disp_x", "Displacement in the x-direction");
  params.addRequiredCoupledVar("disp_y", "Displacement in the y-direction");
  params.addRequiredParam<unsigned int>("component", "Radial = 0, tangential = 1");
  params.addRequiredParam<Point>("center_point",
                                 "Location of the center point of the cylindrical coordinates");
  return params;
}

CylindricalRankOneAux::CylindricalRankOneAux(const InputParameters & parameters)
 : AuxKernel(parameters),
  _disp_x(coupledValue("disp_x")),
  _disp_y(coupledValue("disp_y")), 
  _component(getParam<unsigned int>("component")), 
  _center_point(getParam<Point>("center_point"))
{
}

Real
CylindricalRankOneAux::computeValue()
{
  Point loc_from_center =  *_current_node - _center_point;

  Real theta = std::atan2(loc_from_center(1), loc_from_center(0));
  
  if (_component == 0){
    return (std::cos(theta) * _disp_x[_qp] + std::sin(theta) * _disp_y[_qp]);
  }
  else if (_component == 1){
    return (-std::sin(theta)* _disp_x[_qp] + std::cos(theta)* _disp_y[_qp]);
  }
  else 
    return 0;  
}
