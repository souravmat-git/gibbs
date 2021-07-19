//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TernaryPhaseData.h"
#include "BilinearInterpolation.h"
#include "MooseUtils.h"

registerMooseObject("gibbsApp", TernaryPhaseData);

template <>
InputParameters
validParams<TernaryPhaseData>()
{
  InputParameters params = validParams<ThermoChemicalProperties>();
  params.addRequiredParam<FileName>("table_name", "Phase data in a table");
  params.addClassDescription("Given any tabulated property for a phase..."
                              "this class returns the interpolated values"); 
  return params;
}

TernaryPhaseData::TernaryPhaseData(const InputParameters & parameters)
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
TernaryPhaseData::~TernaryPhaseData() 
{}

void
TernaryPhaseData::initialSetup()
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
      
    //Variable to store the array of free energy
    std::vector<Real> _fenergy = _table_reader.getData(_col_names[2]);
    
    //Variable to store the array of chemical potential of component A
    std::vector<Real> _chem_pot_A = _table_reader.getData(_col_names[3]);
    
    //variable to store the array of chemical potential of component B
    std::vector<Real> _chem_pot_B = _table_reader.getData(_col_names[4]);
    
    //variable to store the array of chemical potential of component C
    std::vector<Real> _chem_pot_C = _table_reader.getData(_col_names[5]);
    
    //variable to store second derivates with respect to B
    std::vector<Real> _tfactor_B = _table_reader.getData(_col_names[6]);
    
    //variable to store second derivates with respect to BC
    std::vector<Real> _tfactor_BC = _table_reader.getData(_col_names[7]);
    
    //variable to store second derivates with respect to C
    std::vector<Real> _tfactor_C = _table_reader.getData(_col_names[8]);
    
    //Declare the size of the matrix, Check PieceWiseBilinearMaterial.C
    
    _GmMatrix.reshape(_num_xC, _num_xB);
    _chemAMatrix.reshape(_num_xC,_num_xB);
    _chemBMatrix.reshape(_num_xC, _num_xB);
    _chemCMatrix.reshape(_num_xC, _num_xB);
    _tfactorBMatrix.reshape(_num_xC,_num_xB);
    _tfactorBCMatrix.reshape(_num_xC, _num_xB);
    _tfactorCMatrix.reshape(_num_xC, _num_xB);
    
    //In case of MooseError "Reference outside of ColumnMajorMatrix bounds"
    //check the size of the matrix 
    //std::cout << _num_xC << _num_xB;
    //Assign the vector elements to the matrix elements
    for(MooseIndex(_num_xC) j=0; j < _num_xC; j++)
      for(MooseIndex(_num_xB) i=0; i < _num_xB; i++)
      {
        _GmMatrix(i,j) = _fenergy[i + _num_xC*j];
        _chemAMatrix(i,j) = _chem_pot_A[i+_num_xC*j];
        _chemBMatrix(i,j) = _chem_pot_B[i + _num_xC*j];
        _chemCMatrix(i,j) = _chem_pot_C[i + _num_xC*j];
        _tfactorBMatrix(i,j) = _tfactor_B[i + _num_xC*j];
        _tfactorBCMatrix(i,j) =_tfactor_BC[i + _num_xC*j];
        _tfactorCMatrix(i,j) = _tfactor_C[i + _num_xC*j];
      }
    
    //Set the data for interpolating free energy
   _interpolate_Gm = libmesh_make_unique<BilinearInterpolation>(_xB, _xC, _GmMatrix);

    //Set the data for interpolating chemical_pot_A
   _interpolate_chem_pot_A = libmesh_make_unique<BilinearInterpolation>(_xB,_xC, _chemAMatrix);
    
    //Set the data for interpolating chemical_pot_B
   _interpolate_chem_pot_B = libmesh_make_unique<BilinearInterpolation>(_xB,_xC ,_chemBMatrix);
    
     //Set the data for interpolating chemical_pot_C
    _interpolate_chem_pot_C = libmesh_make_unique<BilinearInterpolation>(_xB,_xC ,_chemCMatrix);

    //Set the data for interpolating second derivatives
    _interpolate_therm_factor_B = libmesh_make_unique<BilinearInterpolation>(_xB,_xC, _tfactorBMatrix);    
       
    //Set the data for interpolating second derivatives
    _interpolate_therm_factor_BC = libmesh_make_unique<BilinearInterpolation>(_xB,_xC, _tfactorBCMatrix); 
    
    //Set the data for interpolating second derivatives
    _interpolate_therm_factor_C  = libmesh_make_unique<BilinearInterpolation>(_xB,_xC, _tfactorCMatrix); 
}

Real 
TernaryPhaseData::free_energy(const Real& _xB, const Real& _xC) const
{   
    return (_interpolate_Gm->sample(_xB,_xC));   
}

Real
TernaryPhaseData::B_diff_pot(const Real& _xB, const Real& _xC) const
{
  //Return the diffusion potential of component B = chem_pot_B - chem_pot_A
  return (_interpolate_chem_pot_B->sample(_xB,_xC) - _interpolate_chem_pot_A->sample(_xB,_xC));
}

Real
TernaryPhaseData::C_diff_pot(const Real& _xB, const Real& _xC) const
{
  //return 0.0;
  //Return the diffusion potential of component C = chem_pot_C - chem_pot_A
  return (_interpolate_chem_pot_C->sample(_xB,_xC) - _interpolate_chem_pot_A->sample(_xB,_xC));
}

Real
TernaryPhaseData::thermodynamic_factor_B(const Real& _xB, const Real& _xC) const
{
  //return 2.0;
  return (_interpolate_therm_factor_B->sample(_xB,_xC));
}

Real
TernaryPhaseData::thermodynamic_factor_BC(const Real& _xB, const Real& _xC) const
{
  //return 1.0;
  return (_interpolate_therm_factor_BC->sample(_xB,_xC));
}

Real
TernaryPhaseData::thermodynamic_factor_C(const Real& _xB, const Real& _xC) const
{
  //return 4.0;
  return (_interpolate_therm_factor_C->sample(_xB,_xC));
}

