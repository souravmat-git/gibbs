#This is a three-phase problem
#coupling mechanics with diffusion
#It builds the chemical properties
#based on the tabulated approach.
#And then combines with a PRH scheme
#The three phases are
#alpha phase - FCC_A1
#beta phase- FCC_L12
#gamma phase - BCC_B2
#T = 1000K
#Simulation details:
#Interfacial energy = 0.5 J/m2
#interfacial width = 1.5e-6 m
#Length = 300e-6 m
#characteristic_length = 0.05e-6 m

[GlobalParams]
 displacements = 'ux uy'
[]

[Mesh]
  file = anisotropic_strain_jump_2D_ti.e
[]

[Variables]

 [./ux]
   order = FIRST
   family = LAGRANGE
   initial_from_file_var = ux
   initial_from_file_timestep = 9
 [../]

 [./uy]
   order = FIRST
   family = LAGRANGE
   initial_from_file_var = uy
   initial_from_file_timestep = 9
 [../]

 [./x_Al]
   order = FIRST
   family = LAGRANGE
   initial_from_file_var = x_Al
   initial_from_file_timestep = 9
 [../]

 [./mu_Al]
   order = FIRST
   family= LAGRANGE
   initial_from_file_var = mu_Al
   initial_from_file_timestep = 9
 [../]

 [./phi_alpha]
  order = FIRST
  family = LAGRANGE
  initial_from_file_var = phi_alpha
  initial_from_file_timestep = 9
 [../]

 [./phi_beta]
  order = FIRST
  family = LAGRANGE
  initial_from_file_var = phi_beta
  initial_from_file_timestep = 9
 [../]

 [./phi_gamma]
  order = FIRST
  family = LAGRANGE
  initial_from_file_var = phi_gamma
  initial_from_file_timestep = 9
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

#[ICs]

  #[./func_phi_alpha]
  #  type = FunctionIC
  #  variable = phi_alpha
  #  function = func_phi_alpha
  #[../]

  #[./func_phi_beta]
  #  type = FunctionIC
  #  variable = phi_beta
  #  function = func_phi_beta
  #[../]

  #[./func_phi_gamma]
  #  type = FunctionIC
  #  variable = phi_gamma
  #  function = func_phi_gamma
  #[../]

  #[./ic_xAl]
  # type = FunctionIC
  # variable = x_Al
  # function = func_xAl
  #[../]

  #[./muAl_IC]
  # type = FunctionIC
  # variable = mu_Al
  # function = func_muAl
  #[../]

#[]

[Functions]

  [./func_phi_alpha]
    type = ParsedFunction
    value = 'if(x<=-250.0, 1.0, 0)'
  [../]

  [./func_phi_beta]
    type = ParsedFunction
    value = 'if(x>-250.0&x<=250.0, 1.0, 0)'
  [../]

  [./func_phi_gamma]
    type = ParsedFunction
    value = 'if(x>250, 1.0, 0)'
  [../]

  #Equilibrium mole fraction of Al in
  # FCC_A1 (alpha) =
  # FCC_L12(beta)  =
  # FCC_L12 (gamma) =

  [./func_xAl]
    type = ParsedFunction
    vars = 'alpha beta gamma'
    vals = 'func_phi_alpha  func_phi_beta func_phi_gamma'
    #value = '1.24644E-1*alpha + ((4.06459E-1 - 2.37162E-1)*beta*x)/1200 + 4.06459E-1*gamma'
    value = 'if(x<=-250.0, 1.24644E-1, if(x>-250.0&x<=250.0,(4.06459E-1-2.37162E-1)*(x/6000) + 2.37162E-1,4.06459E-1))'
  [../]

  [./func_muAl]
    type = ParsedFunction
    vars = 'alpha beta gamma'
    vals = 'func_phi_alpha  func_phi_beta func_phi_gamma'
    value = '-1.63E1*alpha + -1.63E1*beta + -1.022E1*gamma'
  [../]

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



[Kernels]

  #Variable that this kernel operates on
  #is the diffusion potential of comp B

  [./PhaseConc_Al]
    type = MultiPhaseConstraintMu
    variable = mu_Al
    xB = x_Al
    phase_alpha = phi_alpha
    phase_beta  = phi_beta
    phase_gamma = phi_gamma
    xB_gamma    = xB_gamma
    inv_B_tf_gamma = inv_B_tf_gamma
    h_gamma = h_gamma
  [../]

 ######################
 # Diffusion Equation #
 # Component B        #
 ######################

 [./Al_balance]
   type = BinaryMultiPhaseMassBalance
   variable = x_Al
   B_diff_pot = mu_Al
   phase_alpha = phi_alpha
   phase_beta = phi_beta
   phase_gamma = phi_gamma
   h_alpha = h_alpha
   h_beta  = h_beta
   h_gamma = h_gamma
   inv_B_tf_gamma = inv_B_tf_gamma
 [../]

 [./xAl_dot]
   type = TimeDerivative
   variable = x_Al
 [../]


 #######################
 # Allen-Cahn equation #
 # for phase-alpha     #
 #######################

 [./BinaryMultiPhaseDrivingForce_alphabeta]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_alpha     #Phase_1
    phase_2 = phi_beta       #Coupled phase_2
    phase_3 = phi_gamma
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_alpha
    A_chem_pot_2 = A_chem_pot_beta
    xB_1 = xB_alpha
    xB_2 = xB_beta
    dh = dhbeta_dphialpha     #Required by the residual
    d2h = d2hbeta_dphialpha2  #Required by the diag jacobian
    d2h_2 = d2hbeta_dphialpha_dphibeta
    d2h_3 = d2hbeta_dphialpha_dphigamma
    mob_name = L_phi
    nd_factor = nd_factor
  [../]

  [./BinaryMultiPhaseDrivingForce_alphagamma]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_alpha     #Phase_1
    phase_2 = phi_gamma      #Coupled phase_2
    phase_3 = phi_beta       #Coupled phase_3
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_alpha
    A_chem_pot_2 = A_chem_pot_gamma
    xB_1 = xB_alpha
    xB_2 = xB_gamma
    dh = dhgamma_dphialpha     #Required by the residual
    d2h = d2hgamma_dphialpha2  #Required by the diag jacobian
    d2h_2 = d2hgamma_dphialpha_dphigamma
    d2h_3 = d2hgamma_dphialpha_dphibeta
    mob_name = L_phi
    nd_factor = nd_factor
  [../]

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

   [./doublewell_phi_alphabeta]
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

  [./BinaryMultiPhaseDrivingForce_betaalpha]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_beta
    phase_2 = phi_alpha
    phase_3 = phi_gamma
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_beta
    A_chem_pot_2 = A_chem_pot_alpha
    xB_1 = xB_beta
    xB_2 = xB_alpha
    dh = dhalpha_dphibeta
    d2h = d2halpha_dphibeta2
    d2h_2 = d2halpha_dphibeta_dphialpha
    d2h_3 = d2halpha_dphibeta_dphigamma
    mob_name = L_phi
    nd_factor = nd_factor
  [../]

  [./BinaryMultiPhaseDrivingForce_betagamma]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_beta
    phase_2 = phi_gamma
    phase_3 = phi_alpha
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_beta
    A_chem_pot_2 = A_chem_pot_gamma
    xB_1 = xB_beta
    xB_2 = xB_gamma
    dh = dhgamma_dphibeta
    d2h = d2hgamma_dphibeta2
    d2h_2 = d2hgamma_dphibeta_dphigamma
    d2h_3 = d2hgamma_dphibeta_dphialpha
    mob_name = L_phi
    nd_factor = nd_factor
  [../]

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

 [./BinaryMultiPhaseDrivingForce_gammabeta]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_gamma
    phase_2 = phi_beta
    phase_3 = phi_alpha
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_gamma
    A_chem_pot_2 = A_chem_pot_beta
    xB_1 = xB_gamma
    xB_2 = xB_beta
    dh = dhbeta_dphigamma
    d2h = d2hbeta_dphigamma2
    d2h_2 = d2hbeta_dphigamma_dphibeta
    d2h_3 = d2hbeta_dphigamma_dphialpha
    mob_name = L_phi
    nd_factor = nd_factor
  [../]

  [./BinaryMultiPhaseDrivingForce_gammaalpha]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_gamma
    phase_2 = phi_alpha
    phase_3 = phi_beta
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_gamma
    A_chem_pot_2 = A_chem_pot_alpha
    xB_1 = xB_gamma
    xB_2 = xB_alpha
    dh = dhalpha_dphigamma
    d2h = d2halpha_dphigamma2
    d2h_2 = d2halpha_dphigamma_dphialpha
    d2h_3 = d2halpha_dphigamma_dphibeta
    mob_name = L_phi
    nd_factor = nd_factor
  [../]

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
   mob_name   = L_phi
 [../]

 [./dphi_gamma_dt]
   type = TimeDerivative
   variable = phi_gamma
 [../]

  #Mechanics kernel
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

 [./ConstantFieldProperties]
   type = GenericConstantMaterial
   prop_names =  'L_phi        kappa    m   nd_factor  ndf_el'
   prop_values = '2.772e-04    112.5    1.0   55.427   27692.308'
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

   #########################
   # alpha  is gamma-FCC   #
   #########################

  [./alpha_elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    base_name = alpha
    youngs_modulus = 2.8528
    poissons_ratio = 0.30
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

  #############################
  # beta phase is gamma_prime #
  #############################

  [./beta_elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    base_name = beta
    youngs_modulus = 2.6
    poissons_ratio = 0.30
    #outputs = exodus
  [../]

  [./beta_eigen]
    type = ComputeEigenstrain
    base_name = beta
    eigenstrain_name = eigen
    eigen_base = '-0.003  0.0   0.0
                   0.0   -0.003 0.0
                   0.0    0.0   0.0'
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

  ###########################
  # Gamma is B2-NiAl phase  #
  ###########################

  [./gamma_elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    base_name = gamma
    shear_modulus = 1.3831
    poissons_ratio = 0.3285
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

 [./BinaryConjugatePhaseMaterial_FCC]
   type = BinaryConjugatePhaseMaterial
   A_chem_pot = A_chem_pot_alpha
   xB = xB_alpha
   inv_B_tf = inv_B_tf_alpha
   B_diff_pot = mu_Al
   molar_volume = 1.0
   table_object = table_FCC_data
   #outputs = exodus
 [../]

 [./BinaryConjugatePhaseMaterial_L12]
   type = BinaryConjugatePhaseMaterial
   A_chem_pot = A_chem_pot_beta
   xB = xB_beta
   inv_B_tf = inv_B_tf_beta
   B_diff_pot = mu_Al
   molar_volume = 1.0
   table_object = table_L12_data
   #outputs = exodus
[../]

 [./BinaryConjugatePhaseMaterial_B2]
   type = BinaryConjugatePhaseMaterial
   A_chem_pot = A_chem_pot_gamma
   xB = xB_gamma
   inv_B_tf = inv_B_tf_gamma
   B_diff_pot = mu_Al
   molar_volume = 1.0
   table_object = table_B2_data
   #outputs = exodus
 [../]

 [./BinaryConjugateKineticMaterial_FCC]
  type = BinaryConjugateKineticMaterial
  L_BB = L_BB_alpha
  dL_BB_muB = dL_BB_muB_alpha
  B_diff_pot = mu_Al
  table_object = table_FCC_mobility_data
   #outputs = exodus
 [../]


[./BinaryConjugateKineticMaterial_L12]
   type = BinaryConjugateKineticMaterial
   L_BB = L_BB_beta
   dL_BB_muB = dL_BB_muB_beta
   B_diff_pot = mu_Al
   table_object = table_L12_mobility_data
   #outputs = exodus
 [../]

 [./BinaryConstantKineticMaterial_B2]
   type = BinaryConstantKineticMaterial
   L_BB = L_BB_gamma
   dL_BB_muB = dL_BB_muB_gamma
   L_BB_val = 0.022553
   #outputs = exodus
 [../]

 [./ZeroFieldProperties]
     type = GenericConstantMaterial
     prop_names =  'L_BB_delta L_BB_epsilon
                    dL_BB_muB_delta dL_BB_muB_epsilon'
     prop_values = '0 0 0 0'
  [../]

  ###########################
  # Interpolation Function  #
  ###########################

   [./QuantInterpolationFunction]
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
    type = BinaryConjugatePhaseData
    table_name = "../../../NiAl_data/NiAl_invert_data_FCC.csv"
  [../]

  [./table_L12_data]
    type = BinaryConjugatePhaseData
    table_name = "../../../NiAl_data/NiAl_invert_data_L12.csv"
  [../]

   [./table_B2_data]
    type = BinaryConjugatePhaseData
    table_name = "../../../NiAl_data/NiAl_invert_data_B2.csv"
  [../]

  [./table_FCC_mobility_data]
    type = BinaryConjugateMobilityData
    table_name = "../../../NiAl_data/NiAl_mobility_data_FCC_invert.csv"
  [../]

  [./table_L12_mobility_data]
    type = BinaryConjugateMobilityData
    table_name = "../../../NiAl_data/NiAl_mobility_data_L12_invert.csv"
  [../]

[]

[Postprocessors]

  [./computation_time]
   type = PerfGraphData
   section_name = "Root"
   data_type = total
   execute_on = timestep_end
  [../]

  [./pos_A1_L12]
    type = FindValueOnLine
    target = 0.5
    v = phi_alpha
    start_point = '-3000 0 0'
    end_point = '3000 0 0'
    tol = 1e-7
    execute_on = 'initial timestep_end'
 [../]

  [./pos_L12_B2]
    type = FindValueOnLine
    target = 0.5
    v = phi_gamma
    start_point = '-3000 0 0'
    end_point = '3000 0 0'
    tol = 1e-7
    execute_on = 'initial timestep_end'
 [../]

 [./thickness]
   type = DifferencePostprocessor
   value1 = pos_L12_B2
   value2 = pos_A1_L12
   execute_on = 'initial timestep_end'
 [../]
  #[alpha_mechanical_strain_00]
  #   type = ElementAverageValue
  #   variable = alpha_mechanical_strain_00
  #   outputs = exodus
  # []

  #[beta_mechanical_strain_00]
  #   type = ElementAverageValue
  #   variable = beta_mechanical_strain_00
  # []

  #[gamma_mechanical_strain_00]
  #   type = ElementAverageValue
  #   variable = gamma_mechanical_strain_00
  # []

[]

[VectorPostprocessors]
  [./nodal]
   type = NodalValueSampler
   sort_by = id
   variable = 'x_Al ux uy'
  [../]

 [./fields]
   type = LineValueSampler
   num_points = 1201
   start_point = '-3000 0 0'
   end_point = '3000 0 0'
   sort_by = id
   variable = 'x_Al ux uy total_strain_xx total_strain_xy total_strain_yy stress_xx stress_xy stress_yy'
 [../]
[]


[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu strumpack'

  #petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type'
  #petsc_options_value = 'asm       ilu          nonzero'

  l_max_its = 30
  dtmin = 1e-4
  dtmax = 1e6
  l_tol = 1e-4
  nl_max_its = 20
  nl_rel_tol = 1.0e-8
  nl_abs_tol = 1.0e-10

  start_time = 87037.5
  end_time = 9.5665e4

  #num_steps = 10

  [./TimeStepper]
   type = SolutionTimeAdaptiveDT
   dt = 100.0
  [../]

[]

[Outputs]

  file_base = anisotropic_strain_jump_2D_ti_rs
  #exodus = true

  [./exodus]
    type = Exodus
    execute_on = 'initial timestep_end final'
    interval = 2
  [../]

  [./csv]
    type = CSV
    execute_on = 'initial final'
    execute_postprocessors_on = 'initial timestep_end final'
  [../]

  #[./console]
  #  type = Console
  #  execute_on = 'final'
  #  execute_postprocessors_on = final
  #[../]

  [pgraph]
   type = PerfGraphOutput
   execute_on = 'final'
   level = 1
   heaviest_branch = false
  []

[]

[Debug]
  #show_material_props = true
  show_var_residual_norms = true
[]
