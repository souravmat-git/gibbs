//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
#include<iostream>
//MOOSE includes
#include "MultiSmoothCircleIC.h"

//Define a class
class FileGenerator;

template<>
InputParameters validParams<FileGenerator>();

class FileGenerator : public MultiSmoothCircleIC
{
  public:
    FileGenerator(const InputParameters & parameters);
    
    std::vector<Point> _coords;
    std::vector<Real>  _particle_radii;
    virtual Real value(const Point & p) override;

  protected:
   
};

