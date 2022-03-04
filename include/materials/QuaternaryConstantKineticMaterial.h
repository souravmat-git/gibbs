//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

class QuaternaryConstantKineticMaterial;

//MOOSE includes
#include "Material.h"

//template <>
//InputParameters validParams<QuaternaryConstantKineticMaterial>();

class QuaternaryConstantKineticMaterial : public Material
{
public:
  QuaternaryConstantKineticMaterial(const InputParameters & parameters);

  static InputParameters validParams();

private:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

  //String variable to hold the name of the Onsager mobility matrix
  std::string  _L_BB_name, _L_CC_name, _L_DD_name,
                _L_BC_name, _L_BD_name, _L_CD_name,
                _dL_BB_muB_name, _dL_CC_muB_name, _dL_DD_muB_name,
                _dL_BC_muB_name, _dL_BD_muB_name, _dL_CD_muB_name,
                _dL_BB_muC_name, _dL_CC_muC_name, _dL_DD_muC_name,
                _dL_BC_muC_name, _dL_BD_muC_name, _dL_CD_muC_name,
                _dL_BB_muD_name, _dL_CC_muD_name, _dL_DD_muD_name,
                _dL_BC_muD_name, _dL_BD_muD_name, _dL_CD_muD_name;

   //Onsager mobility diagonal terms
   MaterialProperty<Real> & _L_BB;
   MaterialProperty<Real> & _L_CC;
   MaterialProperty<Real> & _L_DD;

   //Onsager mobility L_BC, L_BD, L_CD
   MaterialProperty<Real> & _L_BC;
   MaterialProperty<Real> & _L_BD;
   MaterialProperty<Real> & _L_CD;

   //Lmatrix with respect to muB
   MaterialProperty<Real> & _dL_BB_muB;
   MaterialProperty<Real> & _dL_CC_muB;
   MaterialProperty<Real> & _dL_DD_muB;

   MaterialProperty<Real> & _dL_BC_muB;
   MaterialProperty<Real> & _dL_BD_muB;
   MaterialProperty<Real> & _dL_CD_muB;

   //Lmatrix with respect to muC
   MaterialProperty<Real> & _dL_BB_muC;
   MaterialProperty<Real> & _dL_CC_muC;
   MaterialProperty<Real> & _dL_DD_muC;

   MaterialProperty<Real> & _dL_BC_muC;
   MaterialProperty<Real> & _dL_BD_muC;
   MaterialProperty<Real> & _dL_CD_muC;

   //Lmatrix with respect to muD
   MaterialProperty<Real> & _dL_BB_muD;
   MaterialProperty<Real> & _dL_CC_muD;
   MaterialProperty<Real> & _dL_DD_muD;

   MaterialProperty<Real> & _dL_BC_muD;
   MaterialProperty<Real> & _dL_BD_muD;
   MaterialProperty<Real> & _dL_CD_muD;


   //Value of the coefficients of the Onsager mobility matrix
   const Real & _L_BB_val;
   const Real & _L_CC_val;
   const Real & _L_DD_val;
   const Real & _L_BC_val;
   const Real & _L_BD_val;
   const Real & _L_CD_val;
};
