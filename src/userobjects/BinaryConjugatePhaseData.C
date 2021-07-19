//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BinaryConjugatePhaseData.h"
#include "LinearInterpolation.h"
#include "MooseUtils.h"

registerMooseObject("gibbsApp", BinaryConjugatePhaseData);

template <>
InputParameters
validParams<BinaryConjugatePhaseData>()
{
  InputParameters params = validParams<ThermoChemicalProperties>();
  params.addRequiredParam<FileName>("table_name", "Phase data in a table");
  params.addClassDescription("Given any tabulated property for a phase..."
                              "this class returns the interpolated values"); 
  return params;
}

BinaryConjugatePhaseData::BinaryConjugatePhaseData(const InputParameters & parameters)
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
BinaryConjugatePhaseData::~BinaryConjugatePhaseData() 
{}

void
BinaryConjugatePhaseData::initialSetup()
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
     
    //Check the first column is always diffusion potential of comp B
      if (_col_names[0] != "B_diff_pot") 
        mooseError("First column must be diffusion potential of comp B !!"); 
   }       
   
    //Variable to store the array of B_diff_pot
    std::vector<Real> _B_diff_pot = _table_reader.getData(_col_names[0]);
    
    //Variable to store the array of chemical potential of the dep comp A
    std::vector<Real> _A_chem_pot = _table_reader.getData(_col_names[1]);
    
    //Variable to store the array of mole fraction of component B
    std::vector<Real> _xB = _table_reader.getData(_col_names[2]);
    
    //Variable to store the inverse of the thermodynamic factor
    std::vector<Real> _inv_tfactor_B = _table_reader.getData(_col_names[3]);
    
    //Variable to store the inverse of the third derivative
    //std::vector<Real> _inv_tderivative_B = _table_reader.getData(_col_names[4]);

    //Set the data for interpolating chemical potential of the dependent comp
    _interpolate_chem_pot_A = libmesh_make_unique<LinearInterpolation>(_B_diff_pot, _A_chem_pot);

    //Set the data for interpolating mole fraction of comp B 
    _interpolate_xB = libmesh_make_unique<LinearInterpolation>(_B_diff_pot, _xB);
  
    //Set the data for interpolating the inverse of the thermodynamic factor B
    _interpolate_inv_tf_B = libmesh_make_unique<LinearInterpolation>(_B_diff_pot, _inv_tfactor_B);   
    
    //Set the data for interpolating the inverse of the third derivatives
    //_interpolate_inv_td_B = libmesh_make_unique<LinearInterpolation>(_B_diff_pot, _inv_tderivative_B); 
}

Real 
BinaryConjugatePhaseData::A_chem_pot(const Real& _B_diff_pot) const
{   
    return (_interpolate_chem_pot_A->sample(_B_diff_pot));   
}

Real
BinaryConjugatePhaseData::xB(const Real& _B_diff_pot) const
{
  return (_interpolate_xB->sample(_B_diff_pot));
}

Real
BinaryConjugatePhaseData::inv_therm_factor_B(const Real& _B_diff_pot) const
{
  return (_interpolate_inv_tf_B->sample(_B_diff_pot));
}

//Real
//BinaryConjugatePhaseData::inv_third_deriv_B(const Real& _B_diff_pot) const
//{
//  return (_interpolate_inv_td_B->sample(_B_diff_pot));
//}
