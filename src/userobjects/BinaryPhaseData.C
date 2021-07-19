//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BinaryPhaseData.h"
#include "LinearInterpolation.h"
#include "MooseUtils.h"

registerMooseObject("gibbsApp", BinaryPhaseData);

template <>
InputParameters
validParams<BinaryPhaseData>()
{
  InputParameters params = validParams<ThermoChemicalProperties>();
  params.addRequiredParam<FileName>("table_name", "Phase data in a table");
  params.addClassDescription("Given any tabulated property for a phase..."
                              "this class returns the interpolated values"); 
  return params;
}

BinaryPhaseData::BinaryPhaseData(const InputParameters & parameters)
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
BinaryPhaseData::~BinaryPhaseData() 
{}

void
BinaryPhaseData::initialSetup()
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
     
    //Check the first column is always mole_fraction of component B
      if (_col_names[0] != "xB") 
        mooseError("First column must be mole fraction !!"); 
   }       
   
    //Variable to store the array of mole_fraction of component B
    std::vector<Real> _xB = _table_reader.getData(_col_names[0]);
    
    //Variable to store the array of free energy
    std::vector<Real> _fenergy = _table_reader.getData(_col_names[1]);
    
    //Variable to store the array of diffusion potential of component B
    std::vector<Real> _B_diff_pot = _table_reader.getData(_col_names[2]);
    
    //Variable to store second derivates with respect to B
    std::vector<Real> _B_therm_factor = _table_reader.getData(_col_names[3]);

    //Set the data for interpolating free energy
    _interpolate_free_energy = libmesh_make_unique<LinearInterpolation>(_xB,_fenergy);

    //Set the data for interpolating diffusion potential of comp B    
    _interpolate_B_diff_pot = libmesh_make_unique<LinearInterpolation>(_xB, _B_diff_pot);
  
    //Set the data for interpolating the thermodynamic factor of comp B
    _interpolate_B_therm_factor = libmesh_make_unique<LinearInterpolation>(_xB, _B_therm_factor);    
}

Real 
BinaryPhaseData::free_energy(const Real& _xB) const
{   
    return (_interpolate_free_energy->sample(_xB));   
}

Real
BinaryPhaseData::B_diff_pot(const Real& _xB) const
{
  //Return the diffusion potential of component B = chem_pot_B - chem_pot_A
  return (_interpolate_B_diff_pot->sample(_xB));
}

Real
BinaryPhaseData::thermodynamic_factor(const Real& _xB) const
{
  return (_interpolate_B_therm_factor->sample(_xB));
}
