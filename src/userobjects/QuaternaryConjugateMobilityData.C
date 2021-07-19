//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "QuaternaryConjugateMobilityData.h"
#include "TrilinearInterpolation.h"
#include "MooseUtils.h"

registerMooseObject("gibbsApp", QuaternaryConjugateMobilityData);

template <>
InputParameters
validParams<QuaternaryConjugateMobilityData>()
{
  InputParameters params = validParams<ThermoChemicalProperties>();
  params.addRequiredParam<FileName>("table_name", "Phase data in a table");
  params.addClassDescription("Given any tabulated property for a phase..."
                              "this class returns the interpolated values"); 
  return params;
}

QuaternaryConjugateMobilityData::QuaternaryConjugateMobilityData(const InputParameters & parameters)
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
QuaternaryConjugateMobilityData::~QuaternaryConjugateMobilityData() 
{}

void
QuaternaryConjugateMobilityData::initialSetup()
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
        mooseError("First and second column must be diffusion potential !!"); 
   }       
   
     //Vector variable to store the array of diffusion potentials of components B,C &D
     _B_diff_pot = _table_reader.getData(_col_names[0]);     
     _C_diff_pot = _table_reader.getData(_col_names[1]);     
     _D_diff_pot = _table_reader.getData(_col_names[2]);
    
    
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
    
     //Make the data unique sort-> unique-> resize for comp D
     std::sort(_D_diff_pot.begin(), _D_diff_pot.end());
    _it_D = std::unique(_D_diff_pot.begin(), _D_diff_pot.end());
    _D_diff_pot.resize(std::distance(_D_diff_pot.begin(),_it_D));
    
    //Size of C vector which is assumed to #of colns
    _num_xD = _D_diff_pot.size();
      
    //Variable to store the array of diagonal coeffecients of the mobility L
    std::vector<Real> _LBB = _table_reader.getData(_col_names[3]);
    std::vector<Real> _LCC = _table_reader.getData(_col_names[4]);
    std::vector<Real> _LDD = _table_reader.getData(_col_names[5]);
    
    //Variable to store the array of off-diagonal coeffecints of the mobility L
    std::vector<Real> _LBC = _table_reader.getData(_col_names[6]);
    std::vector<Real> _LBD = _table_reader.getData(_col_names[7]);
    std::vector<Real> _LCD = _table_reader.getData(_col_names[8]);
    
    //variable to store the first derivative of each coeff. with respect to B
    std::vector<Real> _dLBB_muB = _table_reader.getData(_col_names[9]);
    std::vector<Real> _dLCC_muB = _table_reader.getData(_col_names[10]);
    std::vector<Real> _dLDD_muB = _table_reader.getData(_col_names[11]);
   
    std::vector<Real> _dLBC_muB = _table_reader.getData(_col_names[12]);
    std::vector<Real> _dLBD_muB = _table_reader.getData(_col_names[13]);
    std::vector<Real> _dLCD_muB = _table_reader.getData(_col_names[14]);
    
    //variable to store the first derivative of each coeff. with respect to C
    std::vector<Real> _dLBB_muC = _table_reader.getData(_col_names[15]);
    std::vector<Real> _dLCC_muC = _table_reader.getData(_col_names[16]);
    std::vector<Real> _dLDD_muC = _table_reader.getData(_col_names[17]);
   
    std::vector<Real> _dLBC_muC = _table_reader.getData(_col_names[18]);
    std::vector<Real> _dLBD_muC = _table_reader.getData(_col_names[19]);
    std::vector<Real> _dLCD_muC = _table_reader.getData(_col_names[20]);
    
    //variable to store the first derivative of each coeff. with respect to D
    std::vector<Real> _dLBB_muD = _table_reader.getData(_col_names[21]);
    std::vector<Real> _dLCC_muD = _table_reader.getData(_col_names[22]);
    std::vector<Real> _dLDD_muD = _table_reader.getData(_col_names[23]);
   
    std::vector<Real> _dLBC_muD = _table_reader.getData(_col_names[24]);
    std::vector<Real> _dLBD_muD = _table_reader.getData(_col_names[25]);
    std::vector<Real> _dLCD_muD = _table_reader.getData(_col_names[26]);
    
    
    //Set the data for interpolating diagonal terms
   _interpolate_LBB = 
        libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot, _C_diff_pot,_D_diff_pot, _LBB);  
   _interpolate_LCC = 
        libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot, _D_diff_pot, _LCC);       
   _interpolate_LDD = 
         libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot, _C_diff_pot, _D_diff_pot, _LDD);
         
    //Set the data for interpolating off-diagonal terms
   _interpolate_LBC = 
        libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot, _C_diff_pot,_D_diff_pot, _LBC);  
   _interpolate_LBD = 
        libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot, _D_diff_pot, _LBD);
   _interpolate_LCD = 
         libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot, _C_diff_pot, _D_diff_pot, _LCD);
           
    
    //Set the data for interpolating for first derivatives with respect to comp B
    _interpolate_dLBB_muB = 
         libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot, _D_diff_pot,_dLBB_muB);
    _interpolate_dLCC_muB = 
         libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot, _D_diff_pot, _dLCC_muB);    
    _interpolate_dLDD_muB = 
          libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot, _D_diff_pot, _dLDD_muB);
    
    _interpolate_dLBC_muB = 
          libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot, _D_diff_pot,_dLBC_muB);
    _interpolate_dLBD_muB = 
          libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot,_D_diff_pot,_dLBD_muB);
    _interpolate_dLCD_muB = 
          libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot,_D_diff_pot,_dLCD_muB);
      
      //Now...with respect to comp C    
    _interpolate_dLBB_muC = 
         libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot,_D_diff_pot,_dLBB_muC);
    _interpolate_dLCC_muC = 
         libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot,_D_diff_pot,_dLCC_muC);    
    _interpolate_dLDD_muC = 
          libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot,_D_diff_pot,_dLDD_muC);
    
    _interpolate_dLBC_muC = 
          libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot, _D_diff_pot,_dLBC_muC);
    _interpolate_dLBD_muC = 
          libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot,_D_diff_pot,_dLBD_muC);
    _interpolate_dLCD_muC = 
          libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot,_D_diff_pot,_dLCD_muC);
          
     //Now...with respect to comp D       
    _interpolate_dLBB_muD = 
         libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot, _D_diff_pot,_dLBB_muD);
    _interpolate_dLCC_muD = 
         libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot, _D_diff_pot,_dLCC_muD);    
    _interpolate_dLDD_muD = 
          libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot, _D_diff_pot,_dLDD_muD);
    
    _interpolate_dLBC_muD = 
          libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot,_D_diff_pot,_dLBC_muD);
    _interpolate_dLBD_muD = 
          libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot,_D_diff_pot,_dLBD_muD);
    _interpolate_dLCD_muD = 
          libmesh_make_unique<TrilinearInterpolation>(_B_diff_pot,_C_diff_pot,_D_diff_pot,_dLCD_muD);
}

Real 
QuaternaryConjugateMobilityData::L_BB(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{   
    return (_interpolate_LBB->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));   
}

Real
QuaternaryConjugateMobilityData::L_CC(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{
  return (_interpolate_LCC->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));
}

Real
QuaternaryConjugateMobilityData::L_DD(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real & _D_diff_pot) const
{
  return (_interpolate_LDD->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));
}

Real 
QuaternaryConjugateMobilityData::L_BC(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{   
    return (_interpolate_LBC->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));   
}

Real
QuaternaryConjugateMobilityData::L_BD(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{
  return (_interpolate_LBD->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));
}

Real
QuaternaryConjugateMobilityData::L_CD(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real & _D_diff_pot) const
{
  return (_interpolate_LCD->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));
}

Real 
QuaternaryConjugateMobilityData::dL_BB_muB(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{   
    return (_interpolate_dLBB_muB->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));   
}

Real
QuaternaryConjugateMobilityData::dL_CC_muB(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{
  return (_interpolate_dLCC_muB->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));
}

Real
QuaternaryConjugateMobilityData::dL_DD_muB(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{
  return (_interpolate_dLDD_muB->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));
}

Real
QuaternaryConjugateMobilityData::dL_BC_muB(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{
  return (_interpolate_dLBC_muB->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));
}

Real 
QuaternaryConjugateMobilityData::dL_BD_muB(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{   
    return (_interpolate_dLBD_muB->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));   
}

Real
QuaternaryConjugateMobilityData::dL_CD_muB(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{
  return (_interpolate_dLCD_muB->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));
}

Real 
QuaternaryConjugateMobilityData::dL_BB_muC(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{   
    return (_interpolate_dLBB_muC->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));   
}

Real
QuaternaryConjugateMobilityData::dL_CC_muC(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{
  return (_interpolate_dLCC_muC->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));
}

Real
QuaternaryConjugateMobilityData::dL_DD_muC(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{
  return (_interpolate_dLDD_muC->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));
}

Real
QuaternaryConjugateMobilityData::dL_BC_muC(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{
  return (_interpolate_dLBC_muC->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));
}

Real 
QuaternaryConjugateMobilityData::dL_BD_muC(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{   
    return (_interpolate_dLBD_muC->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));   
}

Real
QuaternaryConjugateMobilityData::dL_CD_muC(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{
  return (_interpolate_dLCD_muC->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));
}

Real 
QuaternaryConjugateMobilityData::dL_BB_muD(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{   
    return (_interpolate_dLBB_muD->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));   
}

Real
QuaternaryConjugateMobilityData::dL_CC_muD(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{
  return (_interpolate_dLCC_muD->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));
}

Real
QuaternaryConjugateMobilityData::dL_DD_muD(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{
  return (_interpolate_dLDD_muD->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));
}

Real
QuaternaryConjugateMobilityData::dL_BC_muD(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{
  return (_interpolate_dLBC_muD->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));
}

Real 
QuaternaryConjugateMobilityData::dL_BD_muD(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{   
    return (_interpolate_dLBD_muD->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));   
}

Real
QuaternaryConjugateMobilityData::dL_CD_muD(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const
{
  return (_interpolate_dLCD_muD->sample(_B_diff_pot, _C_diff_pot, _D_diff_pot));
}
