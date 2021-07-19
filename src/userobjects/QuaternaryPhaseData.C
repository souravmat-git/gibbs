//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "QuaternaryPhaseData.h"
#include "TrilinearInterpolation.h"
#include "MooseUtils.h"

registerMooseObject("gibbsApp", QuaternaryPhaseData);

template <>
InputParameters
validParams<QuaternaryPhaseData>()
{
  InputParameters params = validParams<ThermoChemicalProperties>();
  params.addRequiredParam<FileName>("table_name", "Phase data in a table");
  params.addClassDescription("Given any tabulated property for a phase..."
                              "this class returns the interpolated values"); 
  return params;
}

QuaternaryPhaseData::QuaternaryPhaseData(const InputParameters & parameters)
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
QuaternaryPhaseData::~QuaternaryPhaseData() 
{}

void
QuaternaryPhaseData::initialSetup()
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
     
    //Check the first column is always mole_fraction of components B & C
      if (_col_names[0] != "x_B" && _col_names[1] != "x_C" && _col_names[2] != "x_D") 
        mooseError("First and second column must be mole fraction !!"); 
   }       
   
    //Vector variable to store the array of mole_fraction of components B,C & D
     _xB = _table_reader.getData(_col_names[0]);
     
     _xC = _table_reader.getData(_col_names[1]);
     
     _xD = _table_reader.getData(_col_names[2]);
         
    //Make the data unique sort-> unique -> resize for comp B
     std::sort(_xB.begin(), _xB.end());
    _it_B = std::unique(_xB.begin(), _xB.end());
    _xB.resize(std::distance(_xB.begin(),_it_B));
    
    //Size of B vector which is assumed to be #of rows
    _num_xB = _xB.size();
    
    //Make the data unique sort-> unique-> resize for comp C
     std::sort(_xC.begin(), _xC.end());
    _it_C = std::unique(_xC.begin(), _xC.end());
    _xC.resize(std::distance(_xC.begin(),_it_C));
    
    //Make the data unique sort-> unique-> resize for comp D
     std::sort(_xD.begin(), _xD.end());
    _it_D = std::unique(_xD.begin(), _xD.end());
    _xD.resize(std::distance(_xD.begin(),_it_D));
    
    //Size of D vector which is assumed to #of colns
    _num_xD = _xD.size();
         
    //Variable to store the array of chemical potential of component A
    std::vector<Real> _chem_pot_A = _table_reader.getData(_col_names[3]);
    
    //Variable to store the array of chemical potential of component B
    std::vector<Real> _chem_pot_B = _table_reader.getData(_col_names[4]);
    
    //Variable to store the array of chemical potential of component C
    std::vector<Real> _chem_pot_C = _table_reader.getData(_col_names[5]);
    
    //Variable to store the array of chemical potential of component D
    std::vector<Real> _chem_pot_D = _table_reader.getData(_col_names[6]);
    
    //variable to store second derivates with respect to B
    std::vector<Real> _tfactor_B = _table_reader.getData(_col_names[7]);
    
    //variable to store second derivates with respect to C
    std::vector<Real> _tfactor_C = _table_reader.getData(_col_names[8]);
    
     //variable to store second derivates with respect to D
    std::vector<Real> _tfactor_D = _table_reader.getData(_col_names[9]);
    
    //variable to store second derivates with respect to BC
    std::vector<Real> _tfactor_BC = _table_reader.getData(_col_names[10]);
    
    //variable to store second derivates with respect to BD
    std::vector<Real> _tfactor_BD = _table_reader.getData(_col_names[11]);
    
    //variable to store second derivates with respect to CD
    std::vector<Real> _tfactor_CD = _table_reader.getData(_col_names[12]);
    
    //Set the data for interpolating chemical_pot_A
   _interpolate_chem_pot_A = 
           libmesh_make_unique<TrilinearInterpolation>(_xB,_xC, _xD,_chem_pot_A);
           
    //Set the data for interpolating chemical_pot_B
   _interpolate_chem_pot_B = 
           libmesh_make_unique<TrilinearInterpolation>(_xB,_xC, _xD,_chem_pot_B);
           
      //Set the data for interpolating chemical_pot_C
   _interpolate_chem_pot_C = 
           libmesh_make_unique<TrilinearInterpolation>(_xB,_xC, _xD,_chem_pot_C);
           
     //Set the data for interpolating chemical_pot_D
   _interpolate_chem_pot_D = 
           libmesh_make_unique<TrilinearInterpolation>(_xB,_xC, _xD,_chem_pot_D);
           
    //Set the data for interpolating second derivatives mu(B).x(B)
    _interpolate_therm_factor_B =
           libmesh_make_unique<TrilinearInterpolation>(_xB,_xC, _xD,_tfactor_B); 
           
    //Set the data for interpolating second derivatives mu(C).x(C)
    _interpolate_therm_factor_C =
           libmesh_make_unique<TrilinearInterpolation>(_xB,_xC, _xD,_tfactor_C); 
    
    //Set the data for interpolating second derivatives mu(D).x(D)
    _interpolate_therm_factor_D =
           libmesh_make_unique<TrilinearInterpolation>(_xB,_xC, _xD,_tfactor_D);         
    
    //Set the data for interpolating mu(B).x(C)
    _interpolate_therm_factor_BC = 
           libmesh_make_unique<TrilinearInterpolation>(_xB,_xC, _xD, _tfactor_BC);
    
    //Set the data for interpolating mu(B).x(D)
      _interpolate_therm_factor_BD = 
           libmesh_make_unique<TrilinearInterpolation>(_xB,_xC, _xD, _tfactor_BD);
       
    //Set the data for interpolating mu(C).x(D)
    _interpolate_therm_factor_CD =
           libmesh_make_unique<TrilinearInterpolation>(_xB,_xC,_xD, _tfactor_CD); 
    
}

Real
QuaternaryPhaseData::B_diff_pot(const Real& _xB, const Real& _xC, const Real& _xD) const
{
  //Return the diffusion potential of component B = \mu_B - \mu_A 
  return (_interpolate_chem_pot_B->sample(_xB,_xC, _xD) - _interpolate_chem_pot_A->sample(_xB,_xC,_xD));
}

Real
QuaternaryPhaseData::C_diff_pot(const Real& _xB, const Real& _xC, const Real& _xD) const
{
  //Return the diffusion potential of component B = \mu_B - \mu_A 
  return (_interpolate_chem_pot_C->sample(_xB,_xC, _xD) - _interpolate_chem_pot_A->sample(_xB,_xC,_xD));
}

Real
QuaternaryPhaseData::D_diff_pot(const Real& _xB, const Real& _xC, const Real& _xD) const
{
  //Return the diffusion potential of component B = \mu_B - \mu_A 
  return (_interpolate_chem_pot_D->sample(_xB,_xC, _xD) - _interpolate_chem_pot_A->sample(_xB,_xC,_xD));
}

Real
QuaternaryPhaseData::thermodynamic_factor_B(const Real& _xB, const Real& _xC, const Real& _xD) const
{
  //Return the diffusion potential of component B = \mu_B - \mu_A 
  return (_interpolate_therm_factor_B->sample(_xB,_xC,_xD));
}

Real
QuaternaryPhaseData::thermodynamic_factor_C(const Real& _xB, const Real& _xC, const Real& _xD) const
{
  //Return the diffusion potential of component B = \mu_B - \mu_A 
  return (_interpolate_therm_factor_C->sample(_xB,_xC,_xD));
}

Real
QuaternaryPhaseData::thermodynamic_factor_D(const Real& _xB, const Real& _xC, const Real& _xD) const
{
  //Return the diffusion potential of component B = \mu_B - \mu_A 
  return (_interpolate_therm_factor_D->sample(_xB,_xC,_xD));
}


Real
QuaternaryPhaseData::thermodynamic_factor_BC(const Real& _xB, const Real& _xC, const Real & _xD) const
{
  return (_interpolate_therm_factor_BC->sample(_xB,_xC,_xD));
}

Real
QuaternaryPhaseData::thermodynamic_factor_BD(const Real& _xB, const Real& _xC, const Real & _xD) const
{
  return (_interpolate_therm_factor_BD->sample(_xB,_xC, _xD));
}

Real
QuaternaryPhaseData::thermodynamic_factor_CD(const Real& _xB, const Real& _xC, const Real & _xD) const
{
  return (_interpolate_therm_factor_CD->sample(_xB,_xC, _xD));
}
