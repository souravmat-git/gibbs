#This test problem combines PRH
#with grand-potential model
#This is applied for a void
#growth in U02
#The chemical parameters are obtained
#from th paper by Greenquist CMS 2020
#U02-Void two-phase
#Void is called beta phase
#U02 is called alpha phase
#xVa is the concentration of vacancies
#interface_width = 0.20e-6 m
#The ratio of length to breadth is 30:1

[Mesh]
  file = q2D_twophase_UO2_void_ti3.e
[]


[Variables]
  [./eta]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = eta
    initial_from_file_timestep = LATEST #33
  [../]

  [./u_x]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = u_x
    initial_from_file_timestep = LATEST #33
  [../]

  [./u_y]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = u_y
    initial_from_file_timestep = LATEST #33
  [../]

  #Vacanices mole fraction
  [./x_Va]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = x_Va
    initial_from_file_timestep = LATEST #33
  [../]

  #Vacancies diffusion potential
  [./mu_Va]
   order = FIRST
   family = LAGRANGE
   initial_from_file_var = mu_Va
   initial_from_file_timestep = LATEST #33
  [../]
[]

#[ICs]

  #[./ic_eta]
  #  type = FunctionIC
  #  variable = eta
  #  function = ic_eta_fun
  #[../]

  #[./ic_func_xVa]
  #  variable = x_Va
  #  type = FunctionIC
  #  function = ic_func_xVa
  #[../]

#[]

#[Functions]
#  [./ic_eta_fun]
#    type = ParsedFunction
#    value = 'if(x<=0,1,if(x>0,0,e))'
#  [../]

  #Since it is assumed that beta is the void phase
  #The vacancy concentration in the void phase is 1.0
  #The vacancy concentration is based on the equilibrium
  #vacancy concentration

# [./ic_func_xVa]
#    type = ParsedFunction
#    vars = 'xVa_beta xVa_alpha'
#    vals = '9.9999E-1 3.4252E-2'
#    value = 'if(x<=0,xVa_beta,if(x>0, xVa_alpha,e))'
#  [../]
#[]

[AuxVariables]

  [./grad_eta_0]
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

  [./alpha_es_00]
   order = CONSTANT
   family = MONOMIAL
  [../]

  [./beta_es_00]
   order = CONSTANT
   family = MONOMIAL
  [../]

  [./alpha_se]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./beta_se]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [mech_driving_force]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./interface_energy]
    order = CONSTANT
    family = MONOMIAL
  [../]

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

  #alpha_total_strain = beta_total_strain

  [./total_strain_xx]
    type = RankTwoAux
    rank_two_tensor = alpha_total_strain
    index_i = 0
    index_j = 0
    variable = total_strain_xx
  [../]

  [./total_strain_xy]
    type = RankTwoAux
    rank_two_tensor = alpha_total_strain
    index_i = 0
    index_j = 1
    variable = total_strain_xy
  [../]

  [./total_strain_yy]
    type = RankTwoAux
    rank_two_tensor = alpha_total_strain
    index_i = 1
    index_j = 1
    variable = total_strain_yy
  [../]

   [./alpha_es_00]
    type = RankTwoAux
    rank_two_tensor = alpha_elastic_strain
    index_i = 0
    index_j = 0
    variable = alpha_es_00
  [../]

   [./beta_es_00]
    type = RankTwoAux
    rank_two_tensor = beta_elastic_strain
    index_i = 0
    index_j = 0
    variable = beta_es_00
  [../]

  [./alpha_se]
    type = MaterialRealAux
    variable = alpha_se
    property = alpha_strain_energy_density
    execute_on = timestep_end
  [../]

   [./beta_se]
    type = MaterialRealAux
    variable = beta_se
    property = beta_strain_energy_density
    execute_on = timestep_end
  [../]

  [./mech_driving_force]
    type = MechanicalDrivingForce
    variable = mech_driving_force
    strain_jump_name = strain_jump
    nd_factor  = nd_factor
    execute_on = 'initial timestep_end'
  [../]

  [./interface_local_energy]
    type = InterfaceEnergy
    variable = interface_energy
    phase_beta = eta
    barrier_height = m
  [../]

  [./grad_eta_0]
    type = VariableGradientComponent
    variable = grad_eta_0
    component = x
    gradient_variable = eta
  [../]

[]

[BCs]
 [./Periodic]
    [./auto_ux]
      variable = 'u_x'
      auto_direction = 'y'
    [../]

    [./auto_uy]
      variable = 'u_y'
      auto_direction = 'y'
    [../]

    [./auto_xVa]
      variable = 'x_Va'
      auto_direction = 'y'
    [../]

    [./auto_muVa]
      variable = 'mu_Va'
      auto_direction= 'y'
    [../]

    [./auto_eta]
     variable = eta
     auto_direction = 'y'
    [../]
  [../]

  [./left_dispx]
    type =  DirichletBC
    variable = u_x
    boundary = left
    value = 0
  [../]

  [./left_dispy]
    type = DirichletBC
    variable = u_y
    boundary = left
    value = 0
  [../]

  [./right_dispx]
    type = DirichletBC
    variable = u_x
    boundary = right
    value = 5
  [../]

  [./right_dispy]
    type = DirichletBC
    variable = u_y
    boundary = right
    value = -5
  [../]

[]


[Kernels]

  #This kernel calculates the displacement variable
  #provided the strain jump and its derivative
  #with respect to overall strain is supplied

   [./QMomentumBalance3D_X]
     type = QMomentumBalance3D
     variable = u_x
     component = 0
     displacements = 'u_x u_y'
     eta = eta
   [../]

   [./QMomentumBalance3D_Y]
     type = QMomentumBalance3D
     variable = u_y
     component = 1
     displacements = 'u_x u_y'
     eta = eta
   [../]

   #############################
   # Auxillary equation      ###
   #############################

    [./PhaseConc]
     type = BinaryPhaseConstraintMu
     variable = mu_Va
     xB       = x_Va
     eta      = eta
   [../]

   ######################
   # Diffusion Equation #
   ######################

   [./ContinuityEqn_Va]
     type = BinaryMassBalance
     variable = x_Va
     B_diff_pot = mu_Va
     eta = eta
   [../]

   [./xVa_dot]
    type = TimeDerivative
    variable = x_Va
   [../]

   ########################
   # Driving Force due to #
   # difference in GP     #
   ########################


   [./DrivingForceKKS]
     type = DrivingForceKKS
     variable = eta
     B_diff_pot = mu_Va
     mob_name = L_phi
     nd_factor = nd_factor
   [../]


   ########################
   # Allen-Cahn Equations #
   ########################

   #This kernel implements the double well
   #The form of the double-well barrier is
   #f(phi) = phi^(2)(1 - phi^(2))
   #The free energy is explicitly coded
   #within the kernel

    [./DoubleWell]
     type = DoubleWell
     variable = eta
     mob_name = L_phi
     barrier_height = m
   [../]

   [./ACInterface]
     type = ACInterface
     variable = eta
     kappa_name = kappa
     mob_name = L_phi
   [../]

   [./detadt]
     type = TimeDerivative
     variable = eta
   [../]

   ########################
   # Driving traction at  #
   # the interface        #
   ########################

    [./QDrivingForce3D]
     type = QDrivingForce3D
     variable = eta
     strain_jump_name = strain_jump
     mob_name = L_phi
     nd_factor = ndf_el
   [../]

   [./AllenCahnElasticEnergyOffDiag]
     type = AllenCahnElasticEnergyOffDiag
     variable = eta
     displacements = 'u_x u_y'
     mob_name = L_phi
   [../]

   #This kernel is not fully coupled
   #to strain tensor, and hence
   #the Jacobian corresponding
   #to strain is missing

   #[./QDrivingForceGradPhi3D]
   #  type = QDrivingForceGradPhi3D
   #  variable = eta
   #  mob_name = L_phi
   #  nd_factor = nd_factor
   #[../]
[]



[Materials]

  [./consts]
    type = GenericConstantMaterial
    prop_names  = 'L_phi        kappa   m    nd_factor  ndf_el'
    prop_values = '1.077679e-4  1800   1.0   2.198759  266.239819'
  [../]

  ###########################
  # Overall stress          #
  ###########################

  [./overall_stress]
    type = ComputeOverallStress
    strain_jump_name = strain_jump
  [../]

  #############################
  # Jacobian for QDrivingForce#
  # and QDrivingForceGradPhi  #
  #############################

  #Jacobian for QDrivingForce
  [./d2Fdcdstrain]
    type = ComputeQDrivingForceOffDiag
    strain_jump_name = strain_jump
    nd_factor = ndf_el
  [../]

  #[./drivingforce_gradphi]
  #  type = DrivingForceGradPhi
  #  eta  = eta
  #  strain_jump_name = strain_jump
  #[../]

  #####################
  #### alpha phase ####
  #####################

  [./alpha_eigen]
    type = ComputeEigenstrain
    base_name   = alpha
    eigenstrain_name = eigen
    eigen_base =    '0.0 0.0 0.0
                     0.0 0.0 0.0
                     0.0 0.0 0.0'
    #outputs = exodus
  [../]

  #Uranium dioxide

  [./alpha_elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    base_name = alpha
    youngs_modulus = 2.6
    poissons_ratio = 0.3
   #outputs = exodus
  [../]

  #The mechanical properties
  #are calculated in the following
  #order: First, the elastic strain
  #is calculated then the
  #elastic stress is calculated

  [./alpha_elastic_strain]
   type = ComputeSmallStrainAlpha
   displacements = 'u_x u_y'
   base_name = alpha
   eigenstrain_names = eigen
   #outputs = all
  [../]

  [./alpha_elastic_stress]
    type = ComputePhaseStress
    base_name = alpha
  [../]

  [./alpha_strain_energy_density]
    type = ElasticEnergyMaterial
    base_name = alpha
    f_name = alpha_strain_energy_density
    args = ''
    #outputs = exodus
  [../]

  #####################
  #### beta phase  ####
  #####################

  [./beta_eigen]
    type = ComputeEigenstrain
    base_name = beta
    eigenstrain_name = eigen
    eigen_base = '0.00  0.0  0.0
                  0.0  0.0   0.0
                  0.0  0.0   0.0'
    #outputs = exodus
  [../]

  ## The void phase has a modulus
  # which is 10e-4  times lower than the UO2 phase

   [./beta_elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    base_name = beta
    youngs_modulus = 2.6e-4
    poissons_ratio = 0.3
    #outputs = exodus
  [../]

  [./beta_elastic_strain]
    type = ComputeSmallStrainBeta
    displacements = 'u_x u_y'
    base_name = beta
    eigenstrain_names = eigen
    #outputs = all
  [../]

  [./beta_stress]
    type = ComputePhaseStress
    base_name = beta
  [../]

  [./beta_strain_energy_density]
    type = ElasticEnergyMaterial
    base_name = beta
    f_name = beta_strain_energy_density
    args = ''
    #outputs = exodus
  [../]

  ############################
  # K, a, [e], da/de, da/dh  #
  ############################

  #These materials are arranged in
  #the following logical order:
  #1). Calculate the inverse K tensor
  #2). Calculate the vector a
  #3). Calculate the [e] tensor
  #4). Calculate the derivatives of vector a

  #Note: All caluculations are limited
  #within the interfacial region

  [./ComputeInverseKtensor]
    type = ComputeInverseKtensor
    inv_K_name = inv_K
    eta = eta
  [../]

  [./AnisotropicStrainJumpMaterialN]
    type = AnisotropicStrainJumpMaterialN
    displacements = 'u_x u_y'
    a_name = a
    eta    = eta
    outputs = exodus
  [../]

  [./PartialRankOneHomogenization]
    type = PartialRankOneHomogenization
    strain_jump_name = strain_jump
    eta = eta
    outputs = all
  [../]

  [./DerivativeStrainJumpVectorN]
   type = DerivativeStrainJumpMaterialN
   ds_de = ds_de
   da_dphi = da_dphi
   eta = eta
   #outputs = all
  [../]

  ###############################################
  #### Requisite properties for diffusion #######
  ###############################################

  #The k for void phase was taken to be 10 times
  #that of the solid phase
  #for the solid phase k_alpha = 2.4452e10
  #which was then divided by ev = 6.1272e8
  #curvature is equal to K
  #and horizontal shift = xVa_eqm

  [./GPParabolicPhaseMaterial_Void]
    type = GPParabolicPhaseMaterial
    A_chem_pot = A_chem_pot_beta
    xB = xB_beta
    inv_B_tf = inv_B_tf_beta
    B_diff_pot = mu_Va
    curvature = 399.0785
    horiz_shift = 9.9999E-1
    #outputs = exodus
  [../]

  [./GPParabolicPhaseMaterial_UO2]
    type = GPParabolicPhaseMaterial
    A_chem_pot = A_chem_pot_alpha
    xB = xB_alpha
    inv_B_tf = inv_B_tf_alpha
    B_diff_pot = mu_Va
    curvature = 39.9079
    horiz_shift = 3.4252E-8
    #outputs = exodus
  [../]

  [./BinaryConstantKineticMaterial_Void]
   type = BinaryConstantKineticMaterial
   L_BB = L_BB_beta
   dL_BB_muB = dL_BB_muB_beta
   L_BB_val = 0.0025
   #outputs = exodus
  [../]

  [./BinaryConstantKineticMaterial_UO2]
    type = BinaryConstantKineticMaterial
    L_BB = L_BB_alpha
    dL_BB_muB = dL_BB_muB_alpha
    L_BB_val = 0.0251
    #outputs = exodus
  [../]

  ###########################
  # Interpolation Function  #
  ###########################

   [./InterpolationFunction]
    type = InterpolationFunction
    eta = eta
   [../]

[]

[Postprocessors]
  #[K_00]
  #  type = ElementAverageValue
  #  variable = K_00
  #[]

  #[K_11]
  #  type = ElementAverageValue
  #  variable = K_11
  #[]

  #[K_22]
  #  type = ElementAverageValue
  #  variable = K_22
  #[]

  #[K_01]
  #  type = ElementAverageValue
  #  variable = K_01
  #[]

  #[K_02]
  #  type = ElementAverageValue
  #  variable = K_02
  #[]

  #[K_12]
  #  type = ElementAverageValue
  #  variable = K_12
  #[]

  #############################
  # Strain Jump Tensor        #
  #############################

  [strain_jump_00]
    type = ElementAverageValue
    variable = strain_jump_00
  []

  #[strain_jump_11]
  #  type = ElementAverageValue
  #  variable = strain_jump_11
  #[]

  #[strain_jump_22]
  #  type = ElementAverageValue
  #  variable = strain_jump_22
  #[]

  [stress_xx]
    type = ElementAverageValue
    variable = stress_xx
  []

  [./computation_time]
    type = PerfGraphData
    section_name = Transient::PicardSolve
    data_type = total
    execute_on = 'timestep_end'
  [../]

  [./position]
    type = FindValueOnLine
    target = 0.5
    v = eta
    start_point = '-4410 0 0'
    end_point = '4410 0 0'
    tol = 1e-8
    execute_on = 'initial timestep_end'
  [../]

  [./integrated_interface_energy]
   type = ElementIntegralVariablePostprocessor
   variable = interface_energy
   execute_on = 'initial timestep_end'
  [../]
[]

[VectorPostprocessors]

   [./nodal]
    type = NodalValueSampler
    sort_by = id
    variable = 'eta x_Va u_x u_y'
    execute_on = 'initial timestep_end final'
   [../]

   [./fields]
    type = LineValueSampler
    num_points = 301
    start_point = '-4410 0 0'
    end_point = '4410 0 0'
    sort_by = id
    variable = 'eta x_Va u_x u_y interface_energy'
   [../]

  [./derived_fields]
    type = LineValueSampler
    num_points = 301
    start_point = '-4410 0 0'
    end_point   = '4410 0 0'
    sort_by = id
    variable = 'total_strain_xx total_strain_xy total_strain_yy stress_xx stress_xy stress_yy'
  [../]

  #[./driving_force]
  #  type = LineValueSampler
  #  num_points = 301
  #  start_point = '-4410 0 0'
  #  end_point   = '4410 0 0'
  #  sort_by = id
  #  variable = 'mech_driving_force'
  #[../]
[]


[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]


[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = 'PJFNK'

  #petsc_options_iname = '-pc_type  -sub_pc_type '
  #petsc_options_value = 'asm       lu'


  #petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type'
  #petsc_options_value = 'asm      ilu          nonzero'

  #petsc_options = '-snes_converged_reason -ksp_converged_reason'
  #petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_mat_solver_package'
  #petsc_options_value = 'lu NONZERO superlu_dist'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu strumpack'

  l_max_its = 30
  dtmin = 1e-4
  dtmax = 1e6
  l_tol = 1e-6
  nl_max_its = 15
  nl_rel_tol = 1.0e-8
  nl_abs_tol = 1.0e-10

  start_time = 3.927e9   #5.49182e6
  end_time   = 4e9       #5.58884715e6

  [./TimeStepper]
   type = SolutionTimeAdaptiveDT
   dt = 10.0
  [../]

  #[./TimeStepper]
    #dt = 1.0
    #type = IterationAdaptiveDT
    #cutback_factor = 0.4
    #growth_factor = 1.5
    #optimal_iterations = 5
  #[../]
[]


[Outputs]

  file_base = q2D_twophase_UO2_void_ti3_rs
  #exodus = true

  [./exodus]
    type = Exodus
    execute_on = 'initial timestep_end final'
    interval = 50
  [../]

  [./csv]
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
