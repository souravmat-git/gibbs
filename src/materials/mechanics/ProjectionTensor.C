//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ProjectionTensor.h"
#include "RankTwoTensor.h"

registerMooseObject("gibbsApp", ProjectionTensor);

//template<>
InputParameters
ProjectionTensor::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Calculates the projection tensor from the PFV");
  params.addCoupledVar("eta", "Phase-field variable");
  params.addParam<MaterialPropertyName>("projection_tensor_name",
                                            "projection_tensor",
                                            "define the projection tensor");
  return params;
}

ProjectionTensor::ProjectionTensor(const InputParameters & parameters)
 : Material(parameters),
  _grad_eta(coupledGradient("eta")),
  _projection_tensor(
      declareProperty<RankTwoTensor>
      (getParam<MaterialPropertyName>("projection_tensor_name")))

{
}

void
ProjectionTensor::computeQpProperties()
{
  auto & _P = _projection_tensor[_qp];
  _P.zero();

  //Gradient components obtained from phase-field
  const Real px = _grad_eta[_qp](0);
  const Real py = _grad_eta[_qp](1);
  const Real pz = _grad_eta[_qp](2);

  //Norm of the gradient vector
  const Real mag = std::sqrt((px*px + py*py + pz*pz));

  //Calculate the unit normal
  const Real nx = px/mag;
  const Real ny = py/mag;
  const Real nz = pz/mag;

   //Diagonal components of the matrix
  _P(0,0) = 1.0 - nx * nx;
  _P(1,1) = 1.0 - ny * ny;
  _P(2,2) = 1.0 - nz * nz;

  //Off-diaonal components of the matrix
  _P(0,1) = nx * ny;
  _P(0,2) = nx * nz;
  _P(1,2) = ny * nz;

  //Due to symmetry
  _P(1,0) = _P(0,1);
  _P(2,0) = _P(0,2);
  _P(2,1) = _P(1,2);

}
