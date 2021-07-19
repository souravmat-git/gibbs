//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "QuaternaryConjugatePhaseData.h"
#include "TrilinearInterpolation.h"
#include "MooseUtils.h"

registerMooseObject("gibbsApp", QuaternaryConjugatePhaseData);

template <>
InputParameters
validParams<QuaternaryConjugatePhaseData>()
{
  InputParameters params = validParams<ThermoChemicalProperties>();
  params.addRequiredParam<FileName>("table_name", "Phase data in a table");
  params.addClassDescription("Given any tabulated property for a phase..."
                              "this class returns the interpolated values"); 
  return params;
}

QuaternaryConjugatePhaseData::QuaternaryConjugatePhaseData(const InputParameters & parameters)
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
QuaternaryConjugatePhaseData::~QuaternaryConjugatePhaseData() 
{}

void
QuaternaryConjugatePhaseData::initialSetup()
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
      if (_col_names[0] != "B_diff_pot" && _col_names[1] != "C_diff_pot") 
        mooseError("First and second column must be diffusion potentials !!"); 
   }       
   
    //Vector variable to store the array of diffusion potentials of components B &C
     _B_diff_pot = _table_reader.getData(_col_names[0]);
     
     _C_diff_pot = _table_reader.getData(_col_names[1]);
     
     _D_diff_pot = _table_reader.getData(_col_names[2]);
    
    //Make the data unique sort-> unique -> resize for comp B
     std::sort(_B_diff_pot.begin(), _B_diff_pot.end());
    _it_B = std::unique(_B_diff_pot.begin(), _B_diff_pot.end());
    _B_diff_pot.resize(std::distance(_B_diff_pot.begin(),_it_B));
    
    //Size of B vector which is assumed to be #of rows
    _num_xB = _B_diff_pot.size();
    
    //Make the data unique sort-> unique-> resize for comp C
     std::sort(_C_diff_pot.begin(), _C_diff_pot.end());
    _it_C = std::unique(_C_diff_pot.begin(), _C_diff_pot.end());
    _C_diff_pot.resize(std::distance(_C_diff_pot.begin(),_it_C));
    
    //Size of C vector which is assumed to #of colns
    _num_xC = _C_diff_pot.size();
    
    //Make the data unique sort-> unique-> resize for comp D
     std::sort(_D_diff_pot.begin(), _D_diff_pot.end());
    _it_D = std::unique(_D_diff_pot.begin(), _D_diff_pot.end());
    _D_diff_pot.resize(std::distance(_D_diff_pot.begin(),_it_D));
    
    //Size of C vector which is assumed to #of colns
    _num_xD = _D_diff_pot.size();
      
    //Variable to store the array of chemical potential of A
    std::vector<Real> _chem_pot_A = _table_reader.getData(_col_names[3]);
    
    //Variable to store the array of mole fraction  of component B
    std::vector<Real> _xB = _table_reader.getData(_col_names[4]);
    //Variable to store the array of mole fraction  of component C
    std::vector<Real> _xC = _table_reader.getData(_col_names[5]);
    //Variable to store the array of mole fraction  of component C
    std::vector<Real> _xD = _table_reader.getData(_col_names[6]);
    
    //variable to store the inverse of second derivates with respect to B
    std::vector<Real> _inv_tfactor_B = _table_reader.getData(_col_names[7]);
    //variable to store the inverse of second derivates with respect to BC
    std::vector<Real> _inv_tfactor_C = _table_reader.getData(_col_names[8]);
    //variable to store the inverse of second derivates with respect to C
    std::vector<Real> _inv_tfactor_D = _table_reader.getData(_col_names[9]);
    //variable to store the inverse of second derivates with respect to BC
    std::vector<Real> _inv_tfactor_BC = _table_reader.getData(_col_names[10]);
    //variable to store the inverse of second derivates with respect to BD
    std::vector<Real> _inv_tfactor_BD = _table_reader.getData(_col_names[11]);
    //variable to store the inverse of second derivates with respect to CD
    std::vector<Real> _inv_tfactor_CD = _table_reader.getData(_col_names[12]);
    
 
    //Set the data for interpolating chemical_pot_A
   _interpolate_chem_pot_A = 
            libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot, _D_diff_pot, _chem_pot_A);
    
    //Set the data for interpolating mole fraction of B
   _interpolate_xB = 
            libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot, _D_diff_pot, _xB);
    
     //Set the data for interpolating mole fraction of C
    _interpolate_xC = 
            libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot, _D_diff_pot, _xC);

    //Set the data for interpolating second derivatives
    _interpolate_xD = 
            libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot, _D_diff_pot, _xD);    
       
    //Set the data for interpolating second derivatives
    _interpolate_inv_tf_B = 
            libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot,_D_diff_pot, _inv_tfactor_B); 
    
    //Set the data for interpolating second derivatives
    _interpolate_inv_tf_C  = 
          libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot, _D_diff_pot, _inv_tfactor_C); 
          
      //Set the data for interpolating second derivatives
    _interpolate_inv_tf_D  = 
          libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot, _D_diff_pot, _inv_tfactor_D); 
          
    _interpolate_inv_tf_BC = 
            libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot,_D_diff_pot, _inv_tfactor_BC); 
    
    //Set the data for interpolating second derivatives
    _interpolate_inv_tf_BD  = 
          libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot, _D_diff_pot, _inv_tfactor_BD); 
          
      //Set the data for interpolating second derivatives
    _interpolate_inv_tf_CD  = 
          libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot, _D_diff_pot, _inv_tfactor_CD); 
}

Real
QuaternaryConjugatePhaseData::A_chem_pot(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{
  //Return the chemical potential of dependent component
  return (_interpolate_chem_pot_A->sample(_B_diff_pot,_C_diff_pot, _D_diff_pot));
}

Real
QuaternaryConjugatePhaseData::xB(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real & _D_diff_pot) const
{
  //return 0.0;
  //Return the mole fraction of compB
  return (_interpolate_xB->sample(_B_diff_pot,_C_diff_pot, _D_diff_pot));
}

Real
QuaternaryConjugatePhaseData::xC(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real & _D_diff_pot) const
{
  //return 0.0;
  //Return the mole fraction of compB
  return (_interpolate_xC->sample(_B_diff_pot,_C_diff_pot, _D_diff_pot));
}

Real
QuaternaryConjugatePhaseData::xD(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real & _D_diff_pot) const
{
  //return 0.0;
  //Return the mole fraction of compB
  return (_interpolate_xD->sample(_B_diff_pot,_C_diff_pot, _D_diff_pot));
}

Real
QuaternaryConjugatePhaseData::inv_therm_factor_B(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real & _D_diff_pot) const
{
  //return the inverse of second derivative;
  return (_interpolate_inv_tf_B->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));
}

Real
QuaternaryConjugatePhaseData::inv_therm_factor_C(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real & _D_diff_pot) const
{
  return (_interpolate_inv_tf_C->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));
}

Real
QuaternaryConjugatePhaseData::inv_therm_factor_D(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{
  return (_interpolate_inv_tf_D->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));
}

Real
QuaternaryConjugatePhaseData::inv_therm_factor_BC(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{
  return (_interpolate_inv_tf_BC->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));
}

Real
QuaternaryConjugatePhaseData::inv_therm_factor_BD(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{
  return (_interpolate_inv_tf_BD->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));
}

Real
QuaternaryConjugatePhaseData::inv_therm_factor_CD(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{
  return (_interpolate_inv_tf_CD->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));
}

