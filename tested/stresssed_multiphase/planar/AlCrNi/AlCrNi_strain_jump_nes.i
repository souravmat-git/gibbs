#This kernel uses the TabulatedConjugatePhaseData
#Material to solve the equations.
#This is a model Al-Cr-Ni alloy problem
#alpha_phase-FCC_A1
#beta_phase -FCC_L12
#gamma_phase-BCC_B2
#T = 1473 K
#Simulation details:
#Interfacial energy = 0.5 J/m2
#interfacial width = 1.0e-6 m
#Length = 200e-6 m
#characteristic_length = (0.2e-6/6.0)m = 0.03333e-6 m
#characteristic_time = 1.021038e-02 sec
#Al_diff_pot_eqm = -8.834881E4 J/mol
#Cr_diff_pot_eqm =  1.410318E4 J/mol
#Nondimensinal length = (Length/characteristic_length)
#No. of finite elements = 6.0*(Length/interface_width)

[GlobalParams]
 displacements = 'ux uy'
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1200
  ny = 10
  xmin = -3000
  xmax =  3000
  ymin = -100
  ymax =  100
  elem_type = QUAD4
[]


[Variables]

[./ux]
  order = FIRST
  family = LAGRANGE
[../]

[./uy]
  order = FIRST
  family = LAGRANGE
[../]

[./phi_alpha]
  order = FIRST
  family = LAGRANGE
[../]

[./phi_beta]
  order = FIRST
  family = LAGRANGE
[../]

[./phi_gamma]
  order = FIRST
  family = LAGRANGE
[../]

[./x_Al]
  order = FIRST
  family = LAGRANGE
[../]

[./mu_Al]
  order = FIRST
  family = LAGRANGE
[../]

[./x_Cr]
  order = FIRST
  family = LAGRANGE
[../]

[./mu_Cr]
  order = FIRST
  family = LAGRANGE
[../]

[]

[AuxVariables]

  [./ax_alphabeta]
    order = CONSTANT
    family = MONOMIAL
  [../]

  #[./nx_betagamma]
  #  order = CONSTANT
  #  family = MONOMIAL
  #[../]

  [./ax_betagamma]
    order = CONSTANT
    family = MONOMIAL
  [../]

  #Note that the variable name
  #cannot be the same as the
  #name of the specific
  #tensor components as defined
  #in MOOSE

  [./total_strain_xx]
   order = CONSTANT
   family = MONOMIAL
  [../]

  [./total_strain_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./total_strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]

  #[./beta_mechanical_strain_xx]
  # order = CONSTANT
  # family = MONOMIAL
  #[../]

  #[./gamma_mechanical_strain_xx]
  # order = CONSTANT
  # family = MONOMIAL
  #[../]

[]

[AuxKernels]

  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
    variable = stress_xx
  [../]

  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    variable = stress_xy
  [../]

  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    variable = stress_yy
  [../]

  [./ax_alphabeta]
      type = MaterialRealVectorValueAux
      variable = ax_alphabeta
      property = avec_alpha_beta
      component = 0
      #execute_on = timestep_end
  [../]

  [./ax_betagamma]
      type = MaterialRealVectorValueAux
      variable = ax_betagamma
      property = avec_beta_gamma
      component = 0
      #execute_on = timestep_end
  [../]

  [./total_strain_xx]
   type = RankTwoAux
   variable = total_strain_xx
   rank_two_tensor = total_strain
   index_i = 0
   index_j = 0
  [../]

  [./total_strain_xy]
    type = RankTwoAux
    rank_two_tensor = total_strain
    index_i = 0
    index_j = 1
    variable = total_strain_xy
  [../]

   [./total_strain_yy]
    type = RankTwoAux
    rank_two_tensor = total_strain
    index_i = 1
    index_j = 1
    variable = total_strain_yy
   [../]

  #[./beta_mechanical_strain_xx]
  # type = RankTwoAux
  # variable = beta_mechanical_strain_xx
  # rank_two_tensor = beta_mechanical_strain
  # index_i = 0
  # index_j = 0
  #[../]

  #[./gamma_mechanical_strain_xx]
  # type = RankTwoAux
  # variable = gamma_mechanical_strain_xx
  # rank_two_tensor = gamma_mechanical_strain
  # index_i = 0
  # index_j = 0
  #[../]
[]

[BCs]

  [./Periodic]
    [./auto]
     variable = 'ux uy x_Al phi_alpha phi_beta phi_gamma mu_Al'
     auto_direction = 'y'
    [../]
  [../]

  [./left_disp_x]
   type = DirichletBC
   variable = 'ux'
   boundary = 'left'
   value = 0.0
  [../]

  [./left_disp_y]
     type = DirichletBC
     variable = 'uy'
     boundary = 'left'
     value = 0.0
  [../]

  [./right_disp_x]
   type = DirichletBC
   variable = 'ux'
   boundary = 'right'
   value = 0.0
  [../]

  [./right_disp_y]
    type = DirichletBC
    variable = 'uy'
    boundary = 'right'
    value = 0.0
  [../]



[]


[ICs]

  [./ic_phi_alpha]
    type = FunctionIC
    variable = phi_alpha
    function = func_phi_alpha
  [../]

  [./ic_phi_beta]
    type = FunctionIC
    variable = phi_beta
    function = func_phi_beta
  [../]

  [./ic_phi_gamma]
    type = FunctionIC
    variable = phi_gamma
    function = func_phi_gamma
  [../]

  [./ic_xAl]
   type = FunctionIC
   variable = x_Al
   function = func_xAl
  [../]

  [./muAl_IC]
   type = FunctionIC
   variable = mu_Al
   function = func_muAl
  [../]

  [./ic_xCr]
   type = FunctionIC
   variable = x_Cr
   function = func_xCr
  [../]

  [./muCr_IC]
   type = FunctionIC
   variable = mu_Cr
   function = func_muCr
  [../]

[]

[Functions]

  [./func_phi_alpha]
    type = ParsedFunction
    value = 'if(x<=-510, 1.0, 0)'
  [../]

  [./func_phi_beta]
    type = ParsedFunction
    value = 'if(x>-510 & x<=510, 1.0, 0)'
  [../]

  [./func_phi_gamma]
    type = ParsedFunction
    value = 'if(x>510, 1.0, 0)'
  [../]

  #Equilibrium mole fraction of Al in
  # FCC_A1 (alpha) = 1.676072E-1
  # FCC_L12(beta)  = 2.208670E-1
  # BCC_B2 (gamma) = 2.911855E-1

  [./func_xAl]
    type = ParsedFunction
    vars = 'alpha beta gamma'
    vals = 'func_phi_alpha  func_phi_beta func_phi_gamma'
    value = '1.56440E-1*alpha + 2.158735E-1*beta + 3.00809E-1*gamma'
  [../]

  #Equilibrium mole fraction of Cr in
  # FCC_A1 (alpha) = 1.466087E-1
  # FCC_L12(beta)  = 7.575249E-2
  # BCC_B2 (gamma) = 6.756164E-2

  [./func_xCr]
    type = ParsedFunction
    vars = 'alpha beta gamma'
    vals = 'func_phi_alpha  func_phi_beta func_phi_gamma'
    value = '1.59810E-1*alpha + 8.126219E-2*beta + 5.23381E-2*gamma'
  [../]

  #In dimensional terms = -8.682992E4 (J/mol)/RT
  [./func_muAl]
    type = ParsedFunction
    vars = 'alpha beta gamma'
    vals = 'func_phi_alpha  func_phi_beta func_phi_gamma'
    value = '-7.38012E0*alpha  -7.281933E0*beta  -7.21648E0*gamma'
  [../]

   #In dimensional terms = 2.446668E4 (J/mol)/RT
  [./func_muCr]
    type = ParsedFunction
    vars = 'alpha beta gamma'
    vals = 'func_phi_alpha  func_phi_beta func_phi_gamma'
    value = '1.16184E0*alpha + 1.161837E0*beta + 9.00014E-1*gamma'
  [../]

[]


[Kernels]

  #Variable that this kernel operates on
  #is the diffusion potential of comp B

  [./PhaseConc_Al]
    type = TCPhaseConstraintMuB
    variable = mu_Al
    xB = x_Al
    C_diff_pot = mu_Cr
    phase_alpha = phi_alpha
    phase_beta = phi_beta
    phase_gamma = phi_gamma
    xB_gamma = xB_gamma
    inv_B_tf_gamma = inv_B_tf_gamma
    inv_BC_tf_gamma = inv_BC_tf_gamma
    h_alpha = h_alpha
    h_beta = h_beta
    h_gamma = h_gamma
  [../]

  [./PhaseConc_Cr]
    type = TCPhaseConstraintMuC
    variable = mu_Cr
    xC = x_Cr
    B_diff_pot = mu_Al
    phase_alpha = phi_alpha
    phase_beta = phi_beta
    phase_gamma = phi_gamma
    xC_gamma = xC_gamma
    inv_C_tf_gamma = inv_C_tf_gamma
    inv_BC_tf_gamma = inv_BC_tf_gamma
    h_alpha = h_alpha
    h_beta = h_beta
    h_gamma = h_gamma
  [../]

  ######################
  # Diffusion Equation #
  # Component B        #
  ######################

  [./Al_balance]
    type = TCContinuityEquationB
    variable = x_Al
    xC = x_Cr
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    phase_alpha = phi_alpha
    phase_beta = phi_beta
    phase_gamma = phi_gamma
    h_alpha = h_alpha
    h_beta = h_beta
    h_gamma = h_gamma
    inv_B_tf_gamma = inv_B_tf_gamma
    inv_BC_tf_gamma = inv_BC_tf_gamma
    inv_C_tf_gamma = inv_C_tf_gamma
    L_BB_gamma = L_BB_gamma
    L_BC_gamma = L_BC_gamma
    L_CC_gamma = L_CC_gamma
    dL_BB_muB_gamma = dL_BB_muB_gamma
    dL_BC_muB_gamma = dL_BC_muB_gamma
    dL_CC_muB_gamma = dL_CC_muB_gamma
    dL_BB_muC_gamma = dL_BB_muC_gamma
    dL_BC_muC_gamma = dL_BC_muC_gamma
    dL_CC_muC_gamma = dL_CC_muC_gamma
  [../]

  [./xAl_dot]
    type = TimeDerivative
    variable = x_Al
  [../]

  ######################
  # Diffusion Equation #
  # Component C        #
  ######################

  [./Cr_balance]
    type = TCContinuityEquationC
    variable = x_Cr
    xB = x_Al
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    phase_alpha = phi_alpha
    phase_beta = phi_beta
    phase_gamma = phi_gamma
    h_alpha = h_alpha
    h_beta = h_beta
    h_gamma = h_gamma
    inv_B_tf_gamma = inv_B_tf_gamma
    inv_BC_tf_gamma = inv_BC_tf_gamma
    inv_C_tf_gamma = inv_C_tf_gamma
    L_BB_gamma = L_BB_gamma
    L_BC_gamma = L_BC_gamma
    L_CC_gamma = L_CC_gamma
    dL_BB_muB_gamma = dL_BB_muB_gamma
    dL_BC_muB_gamma = dL_BC_muB_gamma
    dL_CC_muB_gamma = dL_CC_muB_gamma
    dL_BB_muC_gamma = dL_BB_muC_gamma
    dL_BC_muC_gamma = dL_BC_muC_gamma
    dL_CC_muC_gamma = dL_CC_muC_gamma
  [../]

  [./xCr_dot]
    type = TimeDerivative
    variable = x_Cr
  [../]

  #######################
  # Allen-Cahn equation #
  # for phase-alpha     #
  #######################

  [./MultiCompMultiPhaseDrivingForce_alphabeta]
    type = MultiCompDrivingForce
    variable = phi_alpha     #Phase_1
    phase_2 = phi_beta       #Coupled phase_2
    phase_3 = phi_gamma
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    A_chem_pot_1 = A_chem_pot_alpha
    A_chem_pot_2 = A_chem_pot_beta
    xB_1 = xB_alpha
    xB_2 = xB_beta
    xC_1 = xB_alpha
    xC_2 = xC_beta
    dh = dhbeta_dphialpha     #Required by the residual
    d2h = d2hbeta_dphialpha2  #Required by the diag jacobian
    d2h_2 = d2hbeta_dphialpha_dphibeta
    d2h_3 = d2hbeta_dphialpha_dphigamma
    mob_name = L_phi
    nd_factor = nd_factor
  [../]


  [./MultiCompMultiPhaseDrivingForce_alphagamma]
    type = MultiCompDrivingForce
    variable = phi_alpha     #Phase_1
    phase_2 = phi_gamma      #Coupled phase_2
    phase_3 = phi_beta       #Coupled phase_3
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    A_chem_pot_1 = A_chem_pot_alpha
    A_chem_pot_2 = A_chem_pot_gamma
    xB_1 = xB_alpha
    xB_2 = xB_gamma
    xC_1 = xC_alpha
    xC_2 = xC_gamma
    dh = dhgamma_dphialpha     #Required by the residual
    d2h = d2hgamma_dphialpha2  #Required by the diag jacobian
    d2h_2 = d2hgamma_dphialpha_dphigamma
    d2h_3 = d2hgamma_dphialpha_dphibeta
    mob_name = L_phi
    nd_factor = nd_factor
  [../]

  ##################################
  # Driving force due to mechanical#
  # contributions for phase alpha  #
  ##################################

  [./QDrivingForceAlpha_MP]
   type = QDrivingForceAlpha_MP
   variable = phi_alpha
   phase_beta = phi_beta
   phase_gamma = phi_gamma
   mob_name = L_phi
   nd_factor = ndf_el
  [../]

  [./AllenCahnElasticEnergyOffDiagAlpha]
    type = AllenCahnElasticEnergyOffDiagPhase
    variable = phi_alpha
    base_name = alpha
    mob_name = L_phi
  [../]

  #This kernel implements the double well
  #The form of the double-well barrier is
  #f(phi) = phi^(2)(1 - phi^(2))
  #The free energy is explicitly coded
  #within the kernel

   [./doublewell_phi_alpha]
    type = MultiPhaseDoubleWell
    variable = phi_alpha
    eta2 = phi_beta
    eta3 = phi_gamma
    gamma = 1.5
    mob_name = L_phi
    barrier_height = m
  [../]

  [./Curvature_alpha]
    type = ACInterface
    variable = phi_alpha
    kappa_name = kappa
    mob_name = L_phi
  [../]

  [./dalpha_dt]
    type = TimeDerivative
    variable = phi_alpha
  [../]

  #######################
  # Allen-Cahn equation #
  # for phase-beta      #
  #######################

  [./MultiCompDrivingForce_betaalpha]
    type = MultiCompDrivingForce
    variable = phi_beta
    phase_2 = phi_alpha
    phase_3 = phi_gamma
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    A_chem_pot_1 = A_chem_pot_beta
    A_chem_pot_2 = A_chem_pot_alpha
    xB_1 = xB_beta
    xB_2 = xB_alpha
    xC_1 = xC_beta
    xC_2 = xC_alpha
    dh = dhalpha_dphibeta
    d2h = d2halpha_dphibeta2
    d2h_2 = d2halpha_dphibeta_dphialpha
    d2h_3 = d2halpha_dphibeta_dphigamma
    mob_name = L_phi
    nd_factor = nd_factor
  [../]

  [./MultiCompDrivingForce_betagamma]
    type = MultiCompDrivingForce
    variable = phi_beta
    phase_2 = phi_gamma
    phase_3 = phi_alpha
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    A_chem_pot_1 = A_chem_pot_beta
    A_chem_pot_2 = A_chem_pot_gamma
    xB_1 = xB_beta
    xB_2 = xB_gamma
    xC_1 = xC_beta
    xC_2 = xC_gamma
    dh = dhgamma_dphibeta
    d2h = d2hgamma_dphibeta2
    d2h_2 = d2hgamma_dphibeta_dphigamma
    d2h_3 = d2hgamma_dphibeta_dphialpha
    mob_name = L_phi
    nd_factor = nd_factor
  [../]

  ##################################
  # Driving force due to mechanical#
  # contributions for phase beta   #
  ##################################

  [./QDrivingForceBeta_MP]
   type = QDrivingForceBeta_MP
   variable = phi_beta
   phase_alpha = phi_alpha
   phase_gamma = phi_gamma
   mob_name = L_phi
   nd_factor = ndf_el
  [../]

  [./AllenCahnElasticEnergyOffDiagBeta]
    type = AllenCahnElasticEnergyOffDiagPhase
    variable = phi_beta
    base_name = beta
    mob_name = L_phi
  [../]

  #This kernel implements the double well
  #The form of the double-well barrier is
  #f(phi) = phi^(2)(1 - phi^(2))
  #The free energy is explicitly coded
  #within the kernel

   [./doublewell_phi_beta]
    type = MultiPhaseDoubleWell
    variable = phi_beta
    eta2 = phi_alpha
    eta3 = phi_gamma
    gamma = 1.5
    mob_name = L_phi
    barrier_height = m
  [../]

  [./Curvature_beta]
    type = ACInterface
    variable = phi_beta
    kappa_name = kappa
    mob_name = L_phi
  [../]

  [./dphi_beta_dt]
    type = TimeDerivative
    variable = phi_beta
  [../]

  #######################
  # Allen-Cahn equation #
  # for phase-gamma     #
  #######################

  [./MultiCompDrivingForce_gammabeta]
    type = MultiCompDrivingForce
    variable = phi_gamma
    phase_2 = phi_beta
    phase_3 = phi_alpha
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    A_chem_pot_1 = A_chem_pot_gamma
    A_chem_pot_2 = A_chem_pot_beta
    xB_1 = xB_gamma
    xB_2 = xB_beta
    xC_1 = xC_gamma
    xC_2 = xC_beta
    dh = dhbeta_dphigamma
    d2h = d2hbeta_dphigamma2
    d2h_2 = d2hbeta_dphigamma_dphibeta
    d2h_3 = d2hbeta_dphigamma_dphialpha
    mob_name = L_phi
    nd_factor = nd_factor
  [../]

  [./MultiCompDrivingForce_gammaalpha]
    type = MultiCompDrivingForce
    variable = phi_gamma
    phase_2 = phi_alpha
    phase_3 = phi_beta
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    A_chem_pot_1 = A_chem_pot_gamma
    A_chem_pot_2 = A_chem_pot_alpha
    xB_1 = xB_gamma
    xB_2 = xB_alpha
    xC_1 = xC_gamma
    xC_2 = xC_alpha
    dh = dhalpha_dphigamma
    d2h = d2halpha_dphigamma2
    d2h_2 = d2halpha_dphigamma_dphialpha
    d2h_3 = d2halpha_dphigamma_dphibeta
    mob_name = L_phi
    nd_factor = nd_factor
  [../]

  ##################################
  # Driving force due to mechanical#
  # contributions for phase gamma  #
  ##################################

  [./QDrivingForceGamma_MP]
   type = QDrivingForceGamma_MP
   variable = phi_gamma
   phase_beta = phi_beta
   phase_alpha = phi_alpha
   mob_name = L_phi
   nd_factor = ndf_el
  [../]

  [./AllenCahnElasticEnergyOffDiagGamma]
    type = AllenCahnElasticEnergyOffDiagPhase
    variable = phi_gamma
    base_name = gamma
    mob_name = L_phi
  [../]


  #This kernel implements the double well
  #The form of the double-well barrier is
  #f(phi) = phi^(2)(1 - phi^(2))
  #The free energy is explicitly coded
  #within the kernel

   [./doublewell_phi_gamma]
    type = MultiPhaseDoubleWell
    variable = phi_gamma
    eta2 = phi_beta
    eta3 = phi_alpha
    gamma = 1.5
    mob_name = L_phi
    barrier_height = m
  [../]

  [./Curvature_gamma]
    type = ACInterface
    variable = phi_gamma
    kappa_name = kappa
    mob_name = L_phi
  [../]

  [./dphi_gamma_dt]
    type = TimeDerivative
    variable = phi_gamma
  [../]

  #############################
  # Balance of linear Momentum#
  #along x-direction          #
  #Mechanics kernel
  #############################

  [./Momentum_Balance_inx]
    type = QMomentumBalance3D_MP
    variable = ux
    phase_alpha = phi_alpha
    phase_beta  = phi_beta
    phase_gamma = phi_gamma
    use_displaced_mesh = false
    component = 0
  [../]

  [./Momentum_Balance_iny]
    type = QMomentumBalance3D_MP
    variable = uy
    phase_alpha = phi_alpha
    phase_beta  = phi_beta
    phase_gamma = phi_gamma
    use_displaced_mesh = false
    component = 1
  [../]

[]


[Materials]
  ##Phenomenological coeff###

   [./ConstantFieldProperties]
     type = GenericConstantMaterial
     prop_names =  'L_phi     kappa    m    nd_factor  ndf_el'
     prop_values = '3.695e-2  112.5   1.0    54.429   18461.538'
  [../]

  #Returns the total strain
  #Assuming no global strain
  #This is equal to mechanical strain
  #This is given as input to
  #the PrerequisiteVector material

  [./total_strain]
   type = TotalStrain
   #outputs = exodus
  [../]

  #This calculates the
  #overall stress and its jacobian

  [./overall_stress]
    type = ComputeOverall3PStress
  [../]

  #The off-diagonal jacobian terms corresponding
  #to alpha, beta and gamma phases

  [./ComputeQDrivingForceOffDiagAlpha_3P]
    type = ComputeQDrivingForceOffDiagAlpha_3P
    nd_factor = ndf_el
  [../]

  [./ComputeQDrivingForceOffDiagBeta_3P]
    type = ComputeQDrivingForceOffDiagBeta_3P
    nd_factor = ndf_el
  [../]

  [./ComputeQDrivingForceOffDiagGamma_3P]
    type = ComputeQDrivingForceOffDiagGamma_3P
    nd_factor = ndf_el
  [../]

  ######################################
  #Elastic properties of each phase is #
  #defined: first the stiffness,       #
  #then the phase strain and stresses  #
  ######################################

   ##################
   # alpha  phase   #
   ##################

   [./alpha_elasticity_tensor]
     type = ComputeElasticityTensor
     base_name = alpha
     C_ijkl    = '2.2409 1.7080 1.7080
                         2.2409 1.7080
                                2.2409
                                      0.9607
                                            0.9607
                                                   0.9607'
    #outputs = exodus
  [../]

  [./alpha_eigen]
    type = ComputeEigenstrain
    base_name = alpha
    eigenstrain_name = eigen
    eigen_base = '0.0  0.0   0.0
                  0.0  0.0   0.0
                  0.0  0.0   0.0'
    #outputs = exodus
  [../]

  [./alpha_elastic_strain]
   type = ComputeSmallStrain3PAlpha
   base_name = alpha
   eigenstrain_names = eigen
   #outputs = all
  [../]

  [./alpha_stress]
    type = ComputePhaseStress
    base_name  = alpha
  [../]

  [./alpha_strain_energy_density]
    type = ElasticEnergyMaterial
    base_name = alpha
    f_name = alpha_strain_energy_density
    args = ''
    #outputs = exodus
  [../]

  ##############
  # beta phase #
  ##############

  [./beta_elasticity_tensor]
   type = ComputeElasticityTensor
   base_name = beta
   C_ijkl = '2.3128 1.6756 1.6756
                    2.3128 1.6756
                           2.3128
                                  1.0
                                      1.0
                                          1.0'
   #outputs = exodus
 [../]

  [./beta_eigen]
    type = ComputeEigenstrain
    base_name = beta
    eigenstrain_name = eigen
    eigen_base = '-0.00  0.0   0.0
                  0.0    -0.00   0.0
                  0.0      0.0   0.0'
    #outputs = exodus
  [../]

  [./beta_elastic_strain]
   type = ComputeSmallStrain3PBeta
   base_name = beta
   eigenstrain_names = eigen
   #outputs = all
  [../]

  [./beta_stress]
    type = ComputePhaseStress
    base_name  = beta
  [../]

  [./beta_strain_energy_density]
    type = ElasticEnergyMaterial
    base_name = beta
    f_name = beta_strain_energy_density
    args = ''
    #outputs = exodus
  [../]

  #####################################
  # Gamma phase  modulus at T= 1473 K #
  #####################################

  [./gamma_elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    base_name = gamma
    shear_modulus = 0.9115
    poissons_ratio = 0.3387
    #outputs = exodus
  [../]

  [./gamma_eigen]
    type = ComputeEigenstrain
    base_name = gamma
    eigenstrain_name = eigen
    eigen_base = '0.00    0.0   0.0
                  0.0     0.0   0.0
                  0.0     0.0   0.0'
    #outputs = exodus
  [../]

  [./gamma_elastic_strain]
   type = ComputeSmallStrain3PGamma
   base_name = gamma
   eigenstrain_names = eigen
   #outputs = all
 [../]

 [./gamma_stress]
   type = ComputePhaseStress
   base_name  = gamma
 [../]

 [./gamma_strain_energy_density]
   type = ElasticEnergyMaterial
   base_name = gamma
   f_name = gamma_strain_energy_density
   args = ''
   #outputs = exodus
 [../]

 ####################################
 # PRH scheme for three phase system#
 ###################################

 #1. Returns m_alpha, m_beta psi_beta, psi_gamma
 #2. All PrerequisiteTensors
 #3. Calculates the jump vectors
 #4. Then calculates the strain jump
 #5. Finally, the derivatives

 [./PrerequisiteVectors]
   type = PrerequisiteVectors
   phase_alpha = phi_alpha
   phase_beta  = phi_beta
   phase_gamma = phi_gamma
   #outputs = exodus
 [../]

 [./PrerequisiteTensors]
   type = PrerequisiteTensors
   phase_alpha = phi_alpha
   phase_beta  = phi_beta
   phase_gamma = phi_gamma
 [../]

 [./DerivativePrerequisiteTensors]
   type = DerivativePrerequisiteTensors
   phase_alpha = phi_alpha
   phase_beta  = phi_beta
   phase_gamma = phi_gamma
 [../]

 [./AnisotropicStrainJump3P]
   type = AnisotropicStrainJump3P
   phase_alpha = phi_alpha
   phase_beta  = phi_beta
   phase_gamma = phi_gamma
   #outputs = exodus
 [../]

 [./MultiPhaseRankOneHomogenization]
   type = MultiPhaseRankOneHomogenization
   phase_alpha = phi_alpha
   phase_beta  = phi_beta
   phase_gamma = phi_gamma
   #outputs = exodus
 [../]

 [./DerivativeStrainJump]
   type = DerivativeStrainJump_3P
   phase_alpha = phi_alpha
   phase_beta  = phi_beta
   phase_gamma = phi_gamma
 [../]

 ###############################################
 #### Requisite properties for diffusion #######
 ###############################################


  [./TernaryConjugatePhaseMaterial_FCC]
    type = TernaryConjugatePhaseMaterial
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    A_chem_pot = A_chem_pot_alpha
    B_mole_fraction = xB_alpha
    C_mole_fraction = xC_alpha
    inv_B_tf   = inv_B_tf_alpha
    inv_BC_tf  = inv_BC_tf_alpha
    inv_C_tf   = inv_C_tf_alpha
    table_object = table_FCC_data
    #outputs = exodus
  [../]

  [./TernaryConjugatePhaseMaterial_L12]
    type = TernaryConjugatePhaseMaterial
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    A_chem_pot = A_chem_pot_beta
    B_mole_fraction = xB_beta
    C_mole_fraction = xC_beta
    inv_B_tf   = inv_B_tf_beta
    inv_BC_tf  = inv_BC_tf_beta
    inv_C_tf   = inv_C_tf_beta
    table_object = table_L12_data
    #outputs = exodus
 [../]

  [./TernaryConjugatePhaseMaterial_B2]
    type = TernaryConjugatePhaseMaterial
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    A_chem_pot = A_chem_pot_gamma
    B_mole_fraction = xB_gamma
    C_mole_fraction = xC_gamma
    inv_B_tf   = inv_B_tf_gamma
    inv_BC_tf  = inv_BC_tf_gamma
    inv_C_tf   = inv_C_tf_gamma
    table_object = table_B2_data
    #outputs = exodus
  [../]

 [./TernaryConjugateKineticMaterial_FCC]
    type = TernaryConjugateKineticMaterial
    L_BB = L_BB_alpha
    L_BC = L_BC_alpha
    L_CC = L_CC_alpha
    dL_BB_muB = dL_BB_muB_alpha
    dL_BC_muB = dL_BC_muB_alpha
    dL_CC_muB = dL_CC_muB_alpha
    dL_BB_muC = dL_BB_muC_alpha
    dL_BC_muC = dL_BC_muC_alpha
    dL_CC_muC = dL_CC_muC_alpha
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    table_object = table_FCC_mobility_data
    #outputs = exodus
 [../]

 [./TernaryConjugateKineticMaterial_L12]
  type = TernaryConjugateKineticMaterial
  L_BB = L_BB_beta
  L_BC = L_BC_beta
  L_CC = L_CC_beta
  dL_BB_muB = dL_BB_muB_beta
  dL_BC_muB = dL_BC_muB_beta
  dL_CC_muB = dL_CC_muB_beta
  dL_BB_muC = dL_BB_muC_beta
  dL_BC_muC = dL_BC_muC_beta
  dL_CC_muC = dL_CC_muC_beta
  B_diff_pot = mu_Al
  C_diff_pot = mu_Cr
  table_object = table_L12_mobility_data
  #outputs = exodus
 [../]

 [./TernaryConstantKineticMaterial]
  type = TernaryConstantKineticMaterial
  L_BB = L_BB_gamma
  L_BC = L_BC_gamma
  L_CC = L_CC_gamma
  dL_BB_muB = dL_BB_muB_gamma
  dL_BC_muB = dL_BC_muB_gamma
  dL_CC_muB = dL_CC_muB_gamma
  dL_BB_muC = dL_BB_muC_gamma
  dL_BC_muC = dL_BC_muC_gamma
  dL_CC_muC = dL_CC_muC_gamma
  L_BB_val = 0.9272
  L_BC_val = 0.0622
  L_CC_val = 0.3022
  #outputs= exodus
 [../]


 [./InterpolationFunction]
   type = QuantInterpolationFunction
   phase_alpha = phi_alpha
   phase_beta  = phi_beta
   phase_gamma = phi_gamma
 [../]

[]

[UserObjects]

  #This userobjcts reads the data from the table
  #for a given diffusion potential
  #and provides the interpolated grand-potential,
  #mole fraction and inverse of the thermodynamic_factor
  #to the material TabulatedPhaseMaterial

  [./table_FCC_data]
    type = TernaryConjugatePhaseData
    table_name = "../../../AlCrNi_data/AlCrNi_invert_data_FCC_A1.csv"
  [../]

 [./table_L12_data]
    type = TernaryConjugatePhaseData
    table_name = ../../../AlCrNi_data/AlCrNi_invert_data_L12.csv
  [../]

  [./table_B2_data]
    type = TernaryConjugatePhaseData
    table_name = ../../../AlCrNi_data/AlCrNi_invert_data_B2.csv
  [../]

  [./table_FCC_mobility_data]
    type = TernaryConjugateMobilityData
    table_name = ../../../AlCrNi_data/AlCrNi_mobility_data_FCC_invert.csv
  [../]

  [./table_L12_mobility_data]
   type = TernaryConjugateMobilityData
   table_name = ../../../AlCrNi_data/AlCrNi_mobility_data_L12_invert.csv
  [../]

[]

[Postprocessors]

  [./computation_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
    execute_on = timestep_end
  [../]

  [./pos_alpha]
    type = FindValueOnLine
    target = 0.5
    v = phi_alpha
    start_point = '-3000 0 0'
    end_point   = ' 3000 0 0'
    tol = 1e-7
    execute_on = 'initial timestep_end final'
  [../]

  [./pos_gamma]
    type = FindValueOnLine
    target = 0.5
    v = phi_gamma
    start_point = '-3000 0 0'
    end_point   = ' 3000 0 0'
    tol = 1e-7
    execute_on = 'initial timestep_end final'
  [../]

  [./Volumefraction_beta]
    type = VolumeFraction
    variable = phi_beta
    execute_on = 'initial timestep_end final'
  [../]

  [./thickness]
    type = DifferencePostprocessor
    value1 = pos_gamma
    value2 = pos_alpha
    execute_on = 'initial timestep_end final'
  [../]

[]

[VectorPostprocessors]

  [./nodal]
   type = NodalValueSampler
   sort_by = id
   variable = 'x_Al x_Cr ux uy'
   execute_on = 'initial final'
  [../]

  [./fields]
    type = LineValueSampler
    num_points = 1201
    start_point = '-3000 0 0'
    end_point = '3000 0 0'
    sort_by = id
    variable = 'x_Al x_Cr ux uy total_strain_xx total_strain_xy total_strain_yy stress_xx stress_xy stress_yy'
    execute_on = 'initial final'
  [../]

[]

[Preconditioning]
  [./full]
    type = SMP
    full = true
  [../]
[]

[Executioner]

  type = Transient
  scheme = bdf2
  solve_type = 'PJFNK'

  #petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type'
  #petsc_options_value = 'asm       ilu          nonzero'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu strumpack'


  l_max_its = 30
  #dtmin = 1e-4
  dtmax = 1e6
  l_tol = 1e-4
  nl_max_its = 25
  nl_rel_tol = 1.0e-8
  nl_abs_tol = 1.0e-10

  #1hr = 3.5258e+05

  end_time = 3.5258e7

  #num_steps = 20

  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt = 1.0 #real_dt = 1sec
  [../]

[]

[Outputs]

  file_base = AlCrNi_strain_jump_nes

  [./exodus]
    type = Exodus
    execute_on = 'initial timestep_end final'
    interval = 2
  [../]

  [./CSV_format]
    type = CSV
    execute_on = 'initial final'
    execute_postprocessors_on = 'initial timestep_end final'
  [../]

  [pgraph]
    type = PerfGraphOutput
    execute_on = 'final'
    level = 1
    heaviest_branch = false
  []

[]

[Debug]
  show_var_residual_norms = true
[]
