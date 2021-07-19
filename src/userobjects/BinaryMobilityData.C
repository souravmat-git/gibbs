//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BinaryMobilityData.h"
#include "LinearInterpolation.h"
#include "MooseUtils.h"

registerMooseObject("gibbsApp", BinaryMobilityData);

template <>
InputParameters
validParams<BinaryMobilityData>()
{
  InputParameters params = validParams<ThermoChemicalProperties>();
  params.addRequiredParam<FileName>("table_name", "Phase data in a table");
  params.addClassDescription("Given any tabulated property for a phase..."
                              "this class returns the interpolated values"); 
  return params;
}

BinaryMobilityData::BinaryMobilityData(const InputParameters & parameters)
  : ThermoChemicalProperties(parameters),
   _table_name(getParam<FileName>("table_name")),
   _table_reader(_table_name, &_communicator),
   _col_names(_table_reader.getNames())
{
  //Lines begining with # are comments
  _table_reader.setComment("#");  
  //tab is a delimiter
  _table_reader.setDelimiter(" ");
}

//This is a destructor
BinaryMobilityData::~BinaryMobilityData() 
{}

void
BinaryMobilityData::initialSetup()
{
  // Check to see if _file_name supplied exists. If it does, that data
  // will be used. If it does not exist, data will be generated and then
  // written to _file_name.
  std::ifstream file(_table_name.c_str());
  if (file.good())
  {
    _console << "Reading tabulated properties from " << _table_name << "\n";
    //read the data in the table and convert into double
    _table_reader.read(); 
     
    //Check the first column is always mole fraction of comp B
      if (_col_names[0] != "xB") 
        mooseError("First column must be mole fraction !!"); 
   }       
   
    //Variable to store the array of B mole fraction
    std::vector<Real> _xB = _table_reader.getData(_col_names[0]);
    
    //Variable to store the array of L_BB
    std::vector<Real> _L_BB  = _table_reader.getData(_col_names[1]);
    
    //Variable to store the array of first derivative of LBB
    std::vector<Real> _dL_BB_xB = _table_reader.getData(_col_names[2]);
    
    //Set the data for interpolating the mobility of component B
    _interpolate_L_BB= libmesh_make_unique<LinearInterpolation>(_xB, _L_BB);

    //Set the data for interpolating mole fraction of comp B 
    _interpolate_dL_BB_xB = libmesh_make_unique<LinearInterpolation>(_xB, _dL_BB_xB);     
}

Real 
BinaryMobilityData::L_BB(const Real& _xB) const
{   
    return (_interpolate_L_BB->sample(_xB));   
}

Real
BinaryMobilityData::dL_BB_xB(const Real& _xB) const
{
  return (_interpolate_dL_BB_xB->sample(_xB));
}
