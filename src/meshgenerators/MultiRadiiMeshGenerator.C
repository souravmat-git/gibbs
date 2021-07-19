//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MultiRadiiMeshGenerator.h"
registerMooseObject("gibbsApp", MultiRadiiMeshGenerator);

InputParameters
MultiRadiiMeshGenerator::validParams()
{
  InputParameters params = ConcentricCircleMeshGenerator::validParams();
  params.addRequiredParam<unsigned int>("num_circles", "Enter the number of circles");
  params.addClassDescription("Generate uniform concentric circles");
  return params;
}

MultiRadiiMeshGenerator::MultiRadiiMeshGenerator(const InputParameters & parameters)
  : ConcentricCircleMeshGenerator(parameters),
    _num_circles(getParam<unsigned int>("num_circles"))
{
   //_radii.resize(1,_num_circles);
   for (unsigned int i=1; i<_num_circles; i++){
     _radii.push_back(_radii[0] + i*_radii[0]);
     _rings.push_back(1);
   }   
}


std::unique_ptr<MeshBase>
MultiRadiiMeshGenerator::generate()
{
  return ConcentricCircleMeshGenerator::generate();
}
