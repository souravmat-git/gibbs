//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FileGenerator.h"

registerMooseObject("gibbsApp", FileGenerator);

template <>
InputParameters validParams<FileGenerator>()
{
  InputParameters params = validParams<MultiSmoothCircleIC>();

  return params;
}
FileGenerator::FileGenerator(const InputParameters & parameters)
  : MultiSmoothCircleIC(parameters)
{
}

Real
FileGenerator::value(const Point &p)
{
    //open a file
    std::ofstream myfile;
    myfile.open("outputfile.dat"); 
    
    myfile << "x y radii\n"; 
    
    //Read the coordinates and the particle radii
    _coords = SmoothCircleBaseIC::_centers;
    _particle_radii = SmoothCircleBaseIC::_radii;
    
    for(unsigned int i=0; i<_particle_radii.size(); i++){
           myfile << std::scientific << _coords[i](0) << " "
                  << std::scientific << _coords[i](1) << " "
                  << std::scientific << _particle_radii[i]<< "\n";
    }
      
 return SmoothCircleBaseIC::value(p);
}
