//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//*This code is based on the code
//*from TabulatedFluidProperties.C

//#ifndef QUATERNARYCONJUGATEKINETICMATERIAL_H
//#define QUATERNARYCONJUGATEKINETICMATERIAL_H

#pragma once

// Forward Declarations
class QuaternaryConjugateKineticMaterial;

//MOOSE includes
#include "TabulatedKineticMaterial.h"
#include "QuaternaryConjugateMobilityData.h"

//template <>
//InputParameters validParams<QuaternaryConjugateKineticMaterial>();

class QuaternaryConjugateKineticMaterial : public TabulatedKineticMaterial
{
public:
  QuaternaryConjugateKineticMaterial(const InputParameters & parameters);

  static InputParameters validParams();

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

    //String variable to hold the ceff, of the mobility
    // matrix obtained from the input file
    std::string _L_BB_name, _L_CC_name, _L_DD_name,
                _L_BC_name, _L_BD_name, _L_CD_name,
                _dL_BB_muB_name, _dL_CC_muB_name, _dL_DD_muB_name,
                _dL_BC_muB_name, _dL_BD_muB_name, _dL_CD_muB_name,
                _dL_BB_muC_name, _dL_CC_muC_name, _dL_DD_muC_name,
                _dL_BC_muC_name, _dL_BD_muC_name, _dL_CD_muC_name,
                _dL_BB_muD_name, _dL_CC_muD_name, _dL_DD_muD_name,
                _dL_BC_muD_name, _dL_BD_muD_name, _dL_CD_muD_name;

   //Onsager mobility diagonal terms
   MaterialProperty<Real> & _L_BB_val;
   MaterialProperty<Real> & _L_CC_val;
   MaterialProperty<Real> & _L_DD_val;

   //Onsager mobility L_BC, L_BD, L_CD
   MaterialProperty<Real> & _L_BC_val;
   MaterialProperty<Real> & _L_BD_val;
   MaterialProperty<Real> & _L_CD_val;

   //Lmatrix with respect to muB
   MaterialProperty<Real> & _dL_BB_muB_val;
   MaterialProperty<Real> & _dL_CC_muB_val;
   MaterialProperty<Real> & _dL_DD_muB_val;

   MaterialProperty<Real> & _dL_BC_muB_val;
   MaterialProperty<Real> & _dL_BD_muB_val;
   MaterialProperty<Real> & _dL_CD_muB_val;

   //Lmatrix with respect to muC
   MaterialProperty<Real> & _dL_BB_muC_val;
   MaterialProperty<Real> & _dL_CC_muC_val;
   MaterialProperty<Real> & _dL_DD_muC_val;

   MaterialProperty<Real> & _dL_BC_muC_val;
   MaterialProperty<Real> & _dL_BD_muC_val;
   MaterialProperty<Real> & _dL_CD_muC_val;

   //Lmatrix with respect to muD
   MaterialProperty<Real> & _dL_BB_muD_val;
   MaterialProperty<Real> & _dL_CC_muD_val;
   MaterialProperty<Real> & _dL_DD_muD_val;

   MaterialProperty<Real> & _dL_BC_muD_val;
   MaterialProperty<Real> & _dL_BD_muD_val;
   MaterialProperty<Real> & _dL_CD_muD_val;

    //Independent variable on which this property depends
   const VariableValue & _B_diff_pot;

    //Independent variable on which this property depends
   const VariableValue & _C_diff_pot;

   //Independent Variable on which this property depends
   const VariableValue & _D_diff_pot;

    //TernaryPhaseData is an userObject, which reads the
    //tables for ternary alloys with the properties of the phase
    //and linearly interpolates the value .
   const QuaternaryConjugateMobilityData & _table_object;

};
//#endif // QUATERNARYCONJUGATEKINETICMATERIAL_H
