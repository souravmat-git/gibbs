//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BinaryConjugateMobilityData.h"
#include "LinearInterpolation.h"
#include "MooseUtils.h"

registerMooseObject("gibbsApp", BinaryConjugateMobilityData);

template <>
InputParameters
validParams<BinaryConjugateMobilityData>()
{
  InputParameters params = validParams<ThermoChemicalProperties>();
  params.addRequiredParam<FileName>("table_name", "Phase data in a table");
  params.addClassDescription("Given any tabulated property for a phase..."
                              "this class returns the interpolated values"); 
  return params;
}

BinaryConjugateMobilityData::BinaryConjugateMobilityData(const InputParameters & parameters)
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
BinaryConjugateMobilityData::~BinaryConjugateMobilityData() 
{}

void
BinaryConjugateMobilityData::initialSetup()
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
     
    //Check the first column is always B_diff_pou
      if (_col_names[0] != "B_diff_pot") 
        mooseError("First column must be mole fraction !!"); 
   }       
   
    //Variable to store the array of B_diff_pot
    std::vector<Real> _B_diff_pot = _table_reader.getData(_col_names[0]);
    
    //Variable to store the array of L_BB
    std::vector<Real> _L_BB  = _table_reader.getData(_col_names[1]);
    
    //Variable to store the array of first derivative of LBB
    std::vector<Real> _dL_BB_muB = _table_reader.getData(_col_names[2]);
    
    //Set the data for interpolating the mobility of component B
    _interpolate_L_BB= libmesh_make_unique<LinearInterpolation>(_B_diff_pot, _L_BB);

    //Set the data for interpolating mole fraction of comp B 
    _interpolate_dL_BB_muB = libmesh_make_unique<LinearInterpolation>(_B_diff_pot, _dL_BB_muB);     
}

Real 
BinaryConjugateMobilityData::L_BB(const Real& _B_diff_pot) const
{   
    return (_interpolate_L_BB->sample(_B_diff_pot));   
}

Real
BinaryConjugateMobilityData::dL_BB_muB(const Real& _B_diff_pot) const
{
  return (_interpolate_dL_BB_muB->sample(_B_diff_pot));
}
