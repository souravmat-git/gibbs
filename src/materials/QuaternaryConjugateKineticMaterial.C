//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "QuaternaryConjugateKineticMaterial.h"

//Only inlcude the header file of the class
//and not its dependency
registerMooseObject("gibbsApp", QuaternaryConjugateKineticMaterial);

template <>
InputParameters
validParams<QuaternaryConjugateKineticMaterial>()
{
  InputParameters params = validParams<TabulatedKineticMaterial>();
  params.addRequiredParam<MaterialPropertyName>("L_BB", "Onsager mobility L_BB");
  params.addRequiredParam<MaterialPropertyName>("L_CC", "Onsager mobility L_CC");
  params.addRequiredParam<MaterialPropertyName>("L_DD", "Onsager mobility L_DD");
  params.addRequiredParam<MaterialPropertyName>("L_BC", "Onsager mobility L_BC");
  params.addRequiredParam<MaterialPropertyName>("L_BD", "Onsager mobility L_BD");
  params.addRequiredParam<MaterialPropertyName>("L_CD", "Onsager mobility L_CD");
  //w.r.t B
  params.addParam<MaterialPropertyName>("dL_BB_muB",0.0,"L_BB w.r.t B");
  params.addParam<MaterialPropertyName>("dL_CC_muB",0.0,"L_CC w.r.t B");
  params.addParam<MaterialPropertyName>("dL_DD_muB",0.0,"L_DD w.r.t B");
  params.addParam<MaterialPropertyName>("dL_BC_muB",0.0,"L_BC w.r.t B");
  params.addParam<MaterialPropertyName>("dL_BD_muB",0.0,"L_BD w.r.t B");
  params.addParam<MaterialPropertyName>("dL_CD_muB",0.0,"L_CD w.r.t B");
  //w.r.t C
  params.addParam<MaterialPropertyName>("dL_BB_muC",0.0,"L_BB w.r.t C");
  params.addParam<MaterialPropertyName>("dL_CC_muC",0.0,"L_CC w.r.t C");
  params.addParam<MaterialPropertyName>("dL_DD_muC",0.0,"L_DD w.r.t C");
  params.addParam<MaterialPropertyName>("dL_BC_muC",0.0,"L_BC w.r.t C");
  params.addParam<MaterialPropertyName>("dL_BD_muC",0.0,"L_BD w.r.t C");
  params.addParam<MaterialPropertyName>("dL_CD_muC",0.0,"L_CD w.r.t C");
  //w.r.t D
  params.addParam<MaterialPropertyName>("dL_BB_muD",0.0,"L_BB w.r.t D");
  params.addParam<MaterialPropertyName>("dL_CC_muD",0.0,"L_CC w.r.t D");
  params.addParam<MaterialPropertyName>("dL_DD_muD",0.0,"L_DD w.r.t D");
  params.addParam<MaterialPropertyName>("dL_BC_muD",0.0,"L_BC w.r.t D");
  params.addParam<MaterialPropertyName>("dL_BD_muD",0.0,"L_BD w.r.t D");
  params.addParam<MaterialPropertyName>("dL_CD_muD",0.0,"L_CD w.r.t D");
  //Coupled variables
  params.addCoupledVar("B_diff_pot", "Diffusion potential of component B");
  params.addCoupledVar("C_diff_pot", "Diffusion potential of component C");
  params.addCoupledVar("D_diff_pot", "Diffusion potential of componenet D");
  params.addRequiredParam<UserObjectName>("table_object", "Name of the table with phase properties");
  params.addClassDescription("Given any tabulated property data for a phase..."
                              "this class returns the interpolated values"); 
  return params;
}

QuaternaryConjugateKineticMaterial::QuaternaryConjugateKineticMaterial(const InputParameters & parameters)
  : TabulatedKineticMaterial(parameters),
   _L_BB_name(getParam<MaterialPropertyName>("L_BB")),
   _L_CC_name(getParam<MaterialPropertyName>("L_CC")),
   _L_DD_name(getParam<MaterialPropertyName>("L_DD")),
   _L_BC_name(getParam<MaterialPropertyName>("L_BC")),
   _L_BD_name(getParam<MaterialPropertyName>("L_BD")),
   _L_CD_name(getParam<MaterialPropertyName>("L_CD")),
   //Derivative with respect to muB   
   _dL_BB_muB_name(getParam<MaterialPropertyName>("dL_BB_muB")),
   _dL_CC_muB_name(getParam<MaterialPropertyName>("dL_CC_muB")),
   _dL_DD_muB_name(getParam<MaterialPropertyName>("dL_DD_muB")),
   _dL_BC_muB_name(getParam<MaterialPropertyName>("dL_BC_muB")),
   _dL_BD_muB_name(getParam<MaterialPropertyName>("dL_BD_muB")),
   _dL_CD_muB_name(getParam<MaterialPropertyName>("dL_CD_muB")),
    //Derivative with respect to muC  
   _dL_BB_muC_name(getParam<MaterialPropertyName>("dL_BB_muC")),
   _dL_CC_muC_name(getParam<MaterialPropertyName>("dL_CC_muC")),
   _dL_DD_muC_name(getParam<MaterialPropertyName>("dL_DD_muC")),
   _dL_BC_muC_name(getParam<MaterialPropertyName>("dL_BC_muC")),
   _dL_BD_muC_name(getParam<MaterialPropertyName>("dL_BD_muC")),
   _dL_CD_muC_name(getParam<MaterialPropertyName>("dL_CD_muC")),
   //Derivative with respect to muD  
   _dL_BB_muD_name(getParam<MaterialPropertyName>("dL_BB_muD")),
   _dL_CC_muD_name(getParam<MaterialPropertyName>("dL_CC_muD")),
   _dL_DD_muD_name(getParam<MaterialPropertyName>("dL_DD_muD")),
   _dL_BC_muD_name(getParam<MaterialPropertyName>("dL_BC_muD")),
   _dL_BD_muD_name(getParam<MaterialPropertyName>("dL_BD_muD")),
   _dL_CD_muD_name(getParam<MaterialPropertyName>("dL_CD_muD")),
    //Declare properties 
   _L_BB_val(declareProperty<Real>(_L_BB_name)),
   _L_CC_val(declareProperty<Real>(_L_CC_name)),
   _L_DD_val(declareProperty<Real>(_L_DD_name)),
   _L_BC_val(declareProperty<Real>(_L_BC_name)),
   _L_BD_val(declareProperty<Real>(_L_BD_name)),
   _L_CD_val(declareProperty<Real>(_L_CD_name)),
    //Derivatives wr.t muB
   _dL_BB_muB_val(declareProperty<Real>(_dL_BB_muB_name)),
   _dL_CC_muB_val(declareProperty<Real>(_dL_CC_muB_name)),
   _dL_DD_muB_val(declareProperty<Real>(_dL_DD_muB_name)),
   _dL_BC_muB_val(declareProperty<Real>(_dL_BC_muB_name)),
   _dL_BD_muB_val(declareProperty<Real>(_dL_BD_muB_name)),
   _dL_CD_muB_val(declareProperty<Real>(_dL_CD_muB_name)),
   //Derivatives wrt muC
   _dL_BB_muC_val(declareProperty<Real>(_dL_BB_muC_name)),
   _dL_CC_muC_val(declareProperty<Real>(_dL_CC_muC_name)),
   _dL_DD_muC_val(declareProperty<Real>(_dL_DD_muC_name)),
   _dL_BC_muC_val(declareProperty<Real>(_dL_BC_muC_name)),
   _dL_BD_muC_val(declareProperty<Real>(_dL_BD_muC_name)),
   _dL_CD_muC_val(declareProperty<Real>(_dL_CD_muC_name)),
   //Derivatives wrt muD
   _dL_BB_muD_val(declareProperty<Real>(_dL_BB_muD_name)),
   _dL_CC_muD_val(declareProperty<Real>(_dL_CC_muD_name)),
   _dL_DD_muD_val(declareProperty<Real>(_dL_DD_muD_name)),
   _dL_BC_muD_val(declareProperty<Real>(_dL_BC_muD_name)),
   _dL_BD_muD_val(declareProperty<Real>(_dL_BD_muD_name)),
   _dL_CD_muD_val(declareProperty<Real>(_dL_CD_muD_name)),
   //Each property is coupled to
   _B_diff_pot(coupledValue("B_diff_pot")),
   _C_diff_pot(coupledValue("C_diff_pot")),
   _D_diff_pot(coupledValue("D_diff_pot")),
   _table_object(getUserObject<QuaternaryConjugateMobilityData>("table_object"))
{
}

void 
QuaternaryConjugateKineticMaterial::computeQpProperties()
{
    //Diagonal terms of the matrix
    _L_BB_val[_qp] = _table_object.L_BB(_B_diff_pot[_qp], _C_diff_pot[_qp], _D_diff_pot[_qp]);
    _L_CC_val[_qp] = _table_object.L_CC(_B_diff_pot[_qp], _C_diff_pot[_qp], _D_diff_pot[_qp]);
    _L_DD_val[_qp] = _table_object.L_DD(_B_diff_pot[_qp], _C_diff_pot[_qp], _D_diff_pot[_qp]);
    
    //off-diagonal terms
    _L_BC_val[_qp] = _table_object.L_BC(_B_diff_pot[_qp], _C_diff_pot[_qp], _D_diff_pot[_qp]);
    _L_BD_val[_qp] = _table_object.L_BD(_B_diff_pot[_qp], _C_diff_pot[_qp], _D_diff_pot[_qp]);
    _L_CD_val[_qp] = _table_object.L_CD(_B_diff_pot[_qp], _C_diff_pot[_qp], _D_diff_pot[_qp]);
        
    //Derivative with respect to muB
    _dL_BB_muB_val[_qp] = _table_object.dL_BB_muB(_B_diff_pot[_qp], _C_diff_pot[_qp], _D_diff_pot[_qp]);
    _dL_CC_muB_val[_qp] = _table_object.dL_CC_muB(_B_diff_pot[_qp], _C_diff_pot[_qp], _D_diff_pot[_qp]);   
    _dL_DD_muB_val[_qp] = _table_object.dL_DD_muB(_B_diff_pot[_qp], _C_diff_pot[_qp], _D_diff_pot[_qp]);
    
    _dL_BC_muB_val[_qp] = _table_object.dL_BC_muB(_B_diff_pot[_qp], _C_diff_pot[_qp], _D_diff_pot[_qp]);
    _dL_BD_muB_val[_qp] = _table_object.dL_BD_muB(_B_diff_pot[_qp], _C_diff_pot[_qp], _D_diff_pot[_qp]);   
    _dL_CD_muB_val[_qp] = _table_object.dL_CD_muB(_B_diff_pot[_qp], _C_diff_pot[_qp], _D_diff_pot[_qp]);
        
    //Derivative with respect to muC
    _dL_BB_muC_val[_qp] = _table_object.dL_BB_muC(_B_diff_pot[_qp], _C_diff_pot[_qp], _D_diff_pot[_qp]);
    _dL_CC_muC_val[_qp] = _table_object.dL_CC_muC(_B_diff_pot[_qp], _C_diff_pot[_qp], _D_diff_pot[_qp]);   
    _dL_DD_muC_val[_qp] = _table_object.dL_DD_muC(_B_diff_pot[_qp], _C_diff_pot[_qp], _D_diff_pot[_qp]);
    
    _dL_BC_muC_val[_qp] = _table_object.dL_BC_muC(_B_diff_pot[_qp], _C_diff_pot[_qp], _D_diff_pot[_qp]);
    _dL_BD_muC_val[_qp] = _table_object.dL_BD_muC(_B_diff_pot[_qp], _C_diff_pot[_qp], _D_diff_pot[_qp]);   
    _dL_CD_muC_val[_qp] = _table_object.dL_CD_muC(_B_diff_pot[_qp], _C_diff_pot[_qp], _D_diff_pot[_qp]);
    
    //Derivative with respect to muD
    _dL_BB_muD_val[_qp] = _table_object.dL_BB_muD(_B_diff_pot[_qp], _C_diff_pot[_qp], _D_diff_pot[_qp]);
    _dL_CC_muD_val[_qp] = _table_object.dL_CC_muD(_B_diff_pot[_qp], _C_diff_pot[_qp], _D_diff_pot[_qp]);   
    _dL_DD_muD_val[_qp] = _table_object.dL_DD_muD(_B_diff_pot[_qp], _C_diff_pot[_qp], _D_diff_pot[_qp]);
    
    _dL_BC_muD_val[_qp] = _table_object.dL_BC_muD(_B_diff_pot[_qp], _C_diff_pot[_qp], _D_diff_pot[_qp]);
    _dL_BD_muD_val[_qp] = _table_object.dL_BD_muD(_B_diff_pot[_qp], _C_diff_pot[_qp], _D_diff_pot[_qp]);   
    _dL_CD_muD_val[_qp] = _table_object.dL_CD_muD(_B_diff_pot[_qp], _C_diff_pot[_qp], _D_diff_pot[_qp]);
}
