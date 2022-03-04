//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//* This Material class supplies the interpolation function,
//* and its first and second derivatives

#include "StrainDependentTaylorApproximation.h"
registerMooseObject("gibbsApp", StrainDependentTaylorApproximation);

//template <>
InputParameters
StrainDependentTaylorApproximation::validParams()
{
  InputParameters params = TabulatedPhaseMaterial::validParams();
  params.addRequiredParam<std::string>("base_name","The phase name");
  params.addRequiredParam<Real>("xB_eqm", "Equilibrium mole fraction in alpha phase");
  params.addRequiredParam<Real>("B_tf_eqm", "Eqm therm. factor of comp B");
  params.addRequiredParam<Real>("B_diff_pot_eqm", "Equilibrium diffusion potential of comp B");
  params.addRequiredParam<Real>("A_chem_pot_eqm", "Equilibrium chemical potential of comp A");
  params.addRequiredCoupledVar("B_diff_pot", "diffusion potential");
  return params;
}

StrainDependentTaylorApproximation::StrainDependentTaylorApproximation(const InputParameters & parameters)
  : TabulatedPhaseMaterial(parameters),
   _base_name(getParam<std::string>("base_name")),
   //constant material parameters
   _xB_eqm(getParam<Real>("xB_eqm")),
   _B_tf_eqm(getParam<Real>("B_tf_eqm")),
   _B_diff_pot_eqm(getParam<Real>("B_diff_pot_eqm")),
   _A_chem_pot_eqm(getParam<Real>("A_chem_pot_eqm")),
   // Derivative of eigenstrains
   _stress(getMaterialProperty<RankTwoTensor>(_base_name + "_stress")),
   _deigen_dmu(getMaterialProperty<RankTwoTensor>(_base_name + "_deigen_dmu")),
   _d2eigen_dmu2(getMaterialProperty<RankTwoTensor>(_base_name + "_d2eigen_dmu2")),
   _stiffness(getMaterialProperty<RankFourTensor>(_base_name + "_elasticity_tensor")),
    //material properties to be computed...
   _xB_val(declareProperty<Real>("xB_" + _base_name )),
   _inv_B_tf_val(declareProperty<Real>("inv_B_tf_" + _base_name)),
   _A_chem_pot_val(declareProperty<Real>("A_chem_pot_" + _base_name)),
    // as functions of coupled variable
   _B_diff_pot(coupledValue("B_diff_pot"))
{
}

void
StrainDependentTaylorApproximation::computeQpProperties(){

    //Mole fraction of B : chemical + elastic part
    Real _xB_chem = (_B_diff_pot[_qp] - (_B_diff_pot_eqm/_Ec))/(_B_tf_eqm/_Ec) + _xB_eqm;
    Real _xB_elastic = (_deigen_dmu[_qp]).doubleContraction(_stress[_qp]);

    //inverse of thermodynamic factor: chemical + elastic part
    Real _inv_B_tf_chem = (1.0/(_B_tf_eqm/_Ec));
    Real _inv_B_tf_elastic = (_d2eigen_dmu2[_qp]).doubleContraction(_stress[_qp])
                            -(_deigen_dmu[_qp]).doubleContraction(_stiffness[_qp]*_deigen_dmu[_qp]);

    //phase comp. alpha as functon of diffusion potential and strain
    _xB_val[_qp] = _xB_chem + _xB_elastic;

    //First derivative of phase composition alpha with respect to mu as function of strain
    _inv_B_tf_val[_qp] = _inv_B_tf_chem +  _inv_B_tf_elastic;

     //Grand-potential is here returned as function of diffusion potential only
    _A_chem_pot_val[_qp] =  ((_A_chem_pot_eqm/_Ec) - _xB_eqm*(_B_diff_pot[_qp] - (_B_diff_pot_eqm/_Ec))
                            -(0.5/(_B_tf_eqm/_Ec))*std::pow((_B_diff_pot[_qp] - (_B_diff_pot_eqm/_Ec)),2.0)) ;

}
