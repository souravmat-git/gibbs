//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef THERMOCHEMICALPROPERTIES_H
//#define THERMOCHEMICALPROPERTIES_H

#pragma once
//Forward declaration
class ThermoChemicalProperties;

//MOOSE inlcudes
#include "GeneralUserObject.h"

template <>
InputParameters validParams<ThermoChemicalProperties>();

//This class is the BaseClasses for reading thermo-chemical data
//And returnng the interpolated data;

class ThermoChemicalProperties : public GeneralUserObject
{

public:
  ThermoChemicalProperties(const InputParameters & parameters);
  
  //Destructor is a member function which has the 
  //same name as the class but they do not take any argument
  //nor return anything source: WSatvich
  virtual ~ThermoChemicalProperties(); //Destructor
  
  virtual void initialize() override {};
  virtual void execute() override {};
  virtual void finalize() override {};

};
//#endif //THERMOCHEMICALPROPERTIES_H
