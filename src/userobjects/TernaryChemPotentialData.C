//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TernaryChemPotentialData.h"
#include "BilinearInterpolation.h"
#include "MooseUtils.h"

registerMooseObject("gibbsApp", TernaryChemPotentialData);

template <>
InputParameters
validParams<TernaryChemPotentialData>()
{
  InputParameters params = validParams<ThermoChemicalProperties>();
  params.addRequiredParam<FileName>("table_name", "Phase data in a table");
  params.addClassDescription("Given any tabulated property for a phase..."
                              "this class returns the interpolated values"); 
  return params;
}

TernaryChemPotentialData::TernaryChemPotentialData(const InputParameters & parameters)
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
TernaryChemPotentialData::~TernaryChemPotentialData() 
{}

void
TernaryChemPotentialData::initialSetup()
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
      if (_col_names[0] != "x_B" && _col_names[1] != "x_C") 
        mooseError("First and second column must be mole fraction !!"); 
   }       
   
    //Vector variable to store the array of mole_fraction of components B &C
     _xB = _table_reader.getData(_col_names[0]);
     
     _xC = _table_reader.getData(_col_names[1]);
    
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
    
    //Size of C vector which is assumed to #of colns
    _num_xC = _xC.size();
         
    //Variable to store the array of chemical potential of component A
    std::vector<Real> _chem_pot_A = _table_reader.getData(_col_names[2]);
    
    //variable to store second derivates with respect to AB
    std::vector<Real> _tfactor_AB = _table_reader.getData(_col_names[3]);
    
    //variable to store second derivates with respect to AC
    std::vector<Real> _tfactor_AC = _table_reader.getData(_col_names[4]);
    
    //Declare the size of the matrix, Check PieceWiseBilinearMaterial.C
    
    _chemAMatrix.reshape(_num_xC,_num_xB);
    _tfactorABMatrix.reshape(_num_xC,_num_xB);
    _tfactorACMatrix.reshape(_num_xC, _num_xB);

    
    //In case of MooseError "Reference outside of ColumnMajorMatrix bounds"
    //check the size of the matrix 
    //std::cout << _num_xC << _num_xB;
    //Assign the vector elements to the matrix elements
    for(MooseIndex(_num_xC) j=0; j < _num_xC; j++)
      for(MooseIndex(_num_xB) i=0; i < _num_xB; i++)
      {
        _chemAMatrix(i,j) = _chem_pot_A[i+_num_xC*j];
        _tfactorABMatrix(i,j) = _tfactor_AB[i + _num_xC*j];
        _tfactorACMatrix(i,j) = _tfactor_AC[i + _num_xC*j];
      }
    
    //Set the data for interpolating chemical_pot_A
   _interpolate_chem_pot_A = libmesh_make_unique<BilinearInterpolation>(_xB,_xC, _chemAMatrix);
    
    //Set the data for interpolating second derivatives
    _interpolate_therm_factor_AB = libmesh_make_unique<BilinearInterpolation>(_xB,_xC, _tfactorABMatrix);    
       
    //Set the data for interpolating second derivatives
    _interpolate_therm_factor_AC = libmesh_make_unique<BilinearInterpolation>(_xB,_xC, _tfactorACMatrix); 
    
}

Real
TernaryChemPotentialData::A_chem_pot(const Real& _xB, const Real& _xC) const
{
  //Return the chemical potential of component A 
  return (_interpolate_chem_pot_A->sample(_xB,_xC));
}


Real
TernaryChemPotentialData::thermodynamic_factor_AB(const Real& _xB, const Real& _xC) const
{
  //return 2.0;
  return (_interpolate_therm_factor_AB->sample(_xB,_xC));
}

Real
TernaryChemPotentialData::thermodynamic_factor_AC(const Real& _xB, const Real& _xC) const
{
  //return 1.0;
  return (_interpolate_therm_factor_AC->sample(_xB,_xC));
}
