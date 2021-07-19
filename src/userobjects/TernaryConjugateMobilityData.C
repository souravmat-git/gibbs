//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TernaryConjugateMobilityData.h"
#include "BilinearInterpolation.h"
#include "MooseUtils.h"

registerMooseObject("gibbsApp", TernaryConjugateMobilityData);

template <>
InputParameters
validParams<TernaryConjugateMobilityData>()
{
  InputParameters params = validParams<ThermoChemicalProperties>();
  params.addRequiredParam<FileName>("table_name", "Phase data in a table");
  params.addClassDescription("Given any tabulated property for a phase..."
                              "this class returns the interpolated values"); 
  return params;
}

TernaryConjugateMobilityData::TernaryConjugateMobilityData(const InputParameters & parameters)
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
TernaryConjugateMobilityData::~TernaryConjugateMobilityData() 
{}

void
TernaryConjugateMobilityData::initialSetup()
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
        mooseError("First and second column must be mole fraction !!"); 
   }       
   
     //Vector variable to store the array of diffusion potentials of components B &C
     _B_diff_pot = _table_reader.getData(_col_names[0]);
     
     _C_diff_pot = _table_reader.getData(_col_names[1]);
    
    
    //Make the data unique sort-> unique -> resize for comp B
     std::sort(_B_diff_pot.begin(), _B_diff_pot.end());
    _it_B = std::unique(_B_diff_pot.begin(),_B_diff_pot.end());
    _B_diff_pot.resize(std::distance(_B_diff_pot.begin(),_it_B));
    
    //Size of B vector which is assumed to be #of rows
    _num_xB = _B_diff_pot.size();
    
    //Make the data unique sort-> unique-> resize for comp C
     std::sort(_C_diff_pot.begin(), _C_diff_pot.end());
    _it_C = std::unique(_C_diff_pot.begin(), _C_diff_pot.end());
    _C_diff_pot.resize(std::distance(_C_diff_pot.begin(),_it_C));
    
    //Size of C vector which is assumed to #of colns
    _num_xC = _C_diff_pot.size();
      
    //Variable to store the array of free energy
    std::vector<Real> _LBB = _table_reader.getData(_col_names[2]);
    
    //Variable to store the array of chemical potential of component A
    std::vector<Real> _LBC = _table_reader.getData(_col_names[3]);
    
    //variable to store the array of chemical potential of component B
    std::vector<Real> _LCC = _table_reader.getData(_col_names[4]);
    
    //variable to store the array of chemical potential of component C
    std::vector<Real> _dLBB_muB = _table_reader.getData(_col_names[5]);
    
    //variable to store second derivates with respect to B
    std::vector<Real> _dLBB_muC = _table_reader.getData(_col_names[6]);
    
    //variable to store second derivates with respect to BC
    std::vector<Real> _dLBC_muB = _table_reader.getData(_col_names[7]);
    
    //variable to store second derivates with respect to C
    std::vector<Real> _dLBC_muC = _table_reader.getData(_col_names[8]);
    
    //variable to store second derivates with respect to BC
    std::vector<Real> _dLCC_muB = _table_reader.getData(_col_names[9]);
    
    //variable to store second derivates with respect to C
    std::vector<Real> _dLCC_muC = _table_reader.getData(_col_names[10]);
    
    
    //Declare the size of the matrix, Check PieceWiseBilinearMaterial.C
    
    _LBBMatrix.reshape(_num_xC, _num_xB);
    _LBCMatrix.reshape(_num_xC,_num_xB);
    _LCCMatrix.reshape(_num_xC, _num_xB);
    //derivative of LBB
    _dLBBmuBMatrix.reshape(_num_xC, _num_xB);
    _dLBBmuCMatrix.reshape(_num_xC,_num_xB);
    //derivatives of LBC
    _dLBCmuBMatrix.reshape(_num_xC, _num_xB);
    _dLBCmuCMatrix.reshape(_num_xC,_num_xB);
    //derivatives of LCC
    _dLCCmuBMatrix.reshape(_num_xC, _num_xB);
    _dLCCmuCMatrix.reshape(_num_xC,_num_xB);
    
    //In case of MooseError "Reference outside of ColumnMajorMatrix bounds"
    //check the size of the matrix 
    //std::cout << _num_xC << _num_xB;
    //Assign the vector elements to the matrix elements
    for(unsigned int j=0; j < _num_xC; j++)
      for(unsigned int i=0; i < _num_xB; i++)
      {
        _LBBMatrix(i,j) = _LBB[i + _num_xC*j];
        _LBCMatrix(i,j) = _LBC[i+_num_xC*j];
        _LCCMatrix(i,j) = _LCC[i + _num_xC*j];
        _dLBBmuBMatrix(i,j) = _dLBB_muB[i + _num_xC*j];
        _dLBBmuCMatrix(i,j) = _dLBB_muC[i + _num_xC*j];
        _dLBCmuBMatrix(i,j) = _dLBC_muB[i + _num_xC*j];
        _dLBCmuCMatrix(i,j) = _dLBC_muC[i + _num_xC*j];
        _dLCCmuBMatrix(i,j) = _dLCC_muB[i + _num_xC*j];
        _dLCCmuCMatrix(i,j) = _dLCC_muC[i + _num_xC*j];
      }
    
    //Set the data for interpolating free energy
   _interpolate_LBB = 
        libmesh_make_unique<BilinearInterpolation>(_B_diff_pot, _C_diff_pot, _LBBMatrix);

    //Set the data for interpolating chemical_pot_A
   _interpolate_LBC = 
        libmesh_make_unique<BilinearInterpolation>(_B_diff_pot,_C_diff_pot, _LBCMatrix);
    
    //Set the data for interpolating chemical_pot_B
   _interpolate_LCC = 
        libmesh_make_unique<BilinearInterpolation>(_B_diff_pot,_C_diff_pot, _LCCMatrix);
    
    //Set the data for interpolating for L_BB derivatives with respect to comp B,C
    _interpolate_dLBB_muB = 
         libmesh_make_unique<BilinearInterpolation>(_B_diff_pot,_C_diff_pot,_dLBBmuBMatrix);

    _interpolate_dLBB_muC = 
         libmesh_make_unique<BilinearInterpolation>(_B_diff_pot,_C_diff_pot,_dLBBmuCMatrix);    
     
    //Set the data for interpolating for L_BC derivatives with respect to compB, C
    _interpolate_dLBC_muB = 
          libmesh_make_unique<BilinearInterpolation>(_B_diff_pot,_C_diff_pot,_dLBCmuBMatrix);
    
    _interpolate_dLBC_muC = 
          libmesh_make_unique<BilinearInterpolation>(_B_diff_pot,_C_diff_pot,_dLBCmuCMatrix);
    
    //Set the data for interpolating for L_CC derivatives with respect to comp C
    _interpolate_dLCC_muB = 
          libmesh_make_unique<BilinearInterpolation>(_B_diff_pot,_C_diff_pot,_dLCCmuBMatrix);
    
    _interpolate_dLCC_muC = 
          libmesh_make_unique<BilinearInterpolation>(_B_diff_pot,_C_diff_pot,_dLCCmuCMatrix);
}

Real 
TernaryConjugateMobilityData::L_BB(const Real& _B_diff_pot, const Real& _C_diff_pot) const
{   
    return (_interpolate_LBB->sample(_B_diff_pot, _C_diff_pot));   
}

Real
TernaryConjugateMobilityData::L_BC(const Real& _B_diff_pot, const Real& _C_diff_pot) const
{
  return (_interpolate_LBC->sample(_B_diff_pot, _C_diff_pot));
}

Real
TernaryConjugateMobilityData::L_CC(const Real& _B_diff_pot, const Real& _C_diff_pot) const
{
  return (_interpolate_LCC->sample(_B_diff_pot, _C_diff_pot));
}

Real 
TernaryConjugateMobilityData::dL_BB_muB(const Real& _B_diff_pot, const Real& _C_diff_pot) const
{   
    return (_interpolate_dLBB_muB->sample(_B_diff_pot, _C_diff_pot));   
}

Real
TernaryConjugateMobilityData::dL_BC_muB(const Real& _B_diff_pot, const Real& _C_diff_pot) const
{
  return (_interpolate_dLBC_muB->sample(_B_diff_pot, _C_diff_pot));
}

Real
TernaryConjugateMobilityData::dL_CC_muB(const Real& _B_diff_pot, const Real& _C_diff_pot) const
{
  return (_interpolate_dLCC_muB->sample(_B_diff_pot, _C_diff_pot));
}

Real
TernaryConjugateMobilityData::dL_BB_muC(const Real& _B_diff_pot, const Real& _C_diff_pot) const
{
  return (_interpolate_dLBB_muC->sample(_B_diff_pot, _C_diff_pot));
}

Real 
TernaryConjugateMobilityData::dL_BC_muC(const Real& _B_diff_pot, const Real& _C_diff_pot) const
{   
    return (_interpolate_dLBC_muC->sample(_B_diff_pot, _C_diff_pot));   
}

Real
TernaryConjugateMobilityData::dL_CC_muC(const Real& _B_diff_pot, const Real& _C_diff_pot) const
{
  return (_interpolate_dLCC_muC->sample(_B_diff_pot, _C_diff_pot));
}
