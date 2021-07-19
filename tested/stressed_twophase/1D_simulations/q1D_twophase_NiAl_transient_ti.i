#This test problem combines PRH
#with grand-potential model
#for a planar problem
#T = 1473 K
#Length = 10e-6 m
#Al_diff_pot = -1.35156e5 J/mol
#interface_width = 0.10e-6 m
#interface_energy = 36.2e-3 J/m2

#The number of finite elements in
#x-direction is based on (Length/delta_x)
#where delta_x = interface_width/6.0

#The ratio of length to breadth is set at
#30:1 and periodic bcs are enforced along y

[Mesh]
  type = GeneratedMesh
  dim  = 1
  nx   = 600
  xmin = -4410
  xmax =  4410
  elem_type = EDGE2
[]


[Variables]
  [./eta]
    order = FIRST
    family = LAGRANGE
  [../]

  [./u_x]
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
[]

[ICs]

  [./ic_eta]
    type = FunctionIC
    variable = eta
    function = ic_eta_fun
  [../]

   [./ic_func_xAl]
    variable = x_Al
    type = FunctionIC
    function = ic_func_xAl
  [../]

[]

[Functions]
  [./ic_eta_fun]
    type = ParsedFunction
    value = 'if(x<=0,1,if(x>0,0,e))'
  [../]

  [./ic_func_xAl]
    type = ParsedFunction
    vars = 'xAl_beta xAl_alpha'
    vals = '2.30E-1 1.79E-1'
    value = 'if(x<=0,xAl_beta,if(x>0, xAl_alpha,e))'
  [../]
[]

[AuxVariables]

  [./grad_eta_0]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./total_strain_xx]
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

  #alpha_total_strain = beta_total_strain

  [./total_strain_xx]
    type = RankTwoAux
    rank_two_tensor = alpha_total_strain
    index_i = 0
    index_j = 0
    variable = total_strain_xx
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

 #[./Periodic]
   #[./auto]
   # variable = 'u_x u_y x_Al eta mu_Al'
   #auto_direction = 'y'
   #[../]
 #[../]

 [./left_disp_x]
  type = DirichletBC
  variable = 'u_x'
  boundary = 'left'
  value = 0.0
 [../]

 [./right_disp_x]
  type = DirichletBC
  variable = 'u_x'
  boundary = 'right'
  value = 0.0
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
    displacements = 'u_x'
    eta = eta
  [../]


  #############################
  # Auxillary equation      ###
  #############################

   [./PhaseConc]
    type = BinaryPhaseConstraintMu
    variable = mu_Al
    xB       = x_Al
    eta      = eta
  [../]

  ######################
  # Diffusion Equation #
  ######################

  [./ContinuityEqn_Al]
    type = BinaryMassBalance
    variable = x_Al
    B_diff_pot = mu_Al
    eta = eta
  [../]

  [./xAl_dot]
   type = TimeDerivative
   variable = x_Al
  [../]

  ########################
  # Driving Force due to #
  # difference in GP     #
  ########################


  [./DrivingForceKKS]
    type = DrivingForceKKS
    variable = eta
    B_diff_pot = mu_Al
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
    displacements = 'u_x'
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
    prop_names  = 'L_phi        kappa    m   nd_factor   ndf_el'
    prop_values = '8.636162e-4  450.0  1.0   25.570802  8673.252556'
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

  #Gamma

  [./alpha_elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    base_name = alpha
    youngs_modulus = 2.8528
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
   displacements = 'u_x'
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
    eigen_base = '-0.003  0.0    0.0
                   0.0   -0.003  0.0
                   0.0   0.0     0.0'
    #outputs = exodus
  [../]

  ## The precipitate phase beta is softer than the matrix
  #For this reason E_beta < E_alpha

   [./beta_elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    base_name = beta
    youngs_modulus = 2.6
    poissons_ratio = 0.3
    #outputs = exodus
  [../]

  [./beta_elastic_strain]
    type = ComputeSmallStrainBeta
    displacements = 'u_x'
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
    displacements = 'u_x'
    a_name = a
    eta    = eta
    #outputs = exodus
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

  #Based on taylor approximation
  #these values are determined

  [./GPtaylorapproximation_L12]
    type = GPTaylorApproximation
    A_chem_pot = A_chem_pot_beta
    xB = xB_beta
    inv_B_tf = inv_B_tf_beta
    B_diff_pot = mu_Al
    char_energy = 1.225e4
    xB_eqm = 2.30730E-1
    B_tf_eqm = 2.86035E5
    B_diff_pot_eqm = -1.08531E5
    A_chem_pot_eqm = -8.56263E4
    #outputs = exodus
  [../]

  [./GPtaylorapproximation_fcc]
    type = GPTaylorApproximation
    A_chem_pot = A_chem_pot_alpha
    xB = xB_alpha
    inv_B_tf = inv_B_tf_alpha
    B_diff_pot = mu_Al
    char_energy = 1.225e4
    xB_eqm = 1.83922E-1
    B_tf_eqm = 3.6927E5
    B_diff_pot_eqm = -1.08531E5
    A_chem_pot_eqm = -8.56263E4
    #outputs = exodus
  [../]

  [./BinaryConstantKineticMaterial_L12]
   type = BinaryConstantKineticMaterial
   L_BB = L_BB_beta
   dL_BB_muB = dL_BB_muB_beta
   L_BB_val = 0.0512
   #outputs = exodus
  [../]

  [./BinaryConstantKineticMaterial_fcc]
    type = BinaryConstantKineticMaterial
    L_BB = L_BB_alpha
    dL_BB_muB = dL_BB_muB_alpha
    L_BB_val = 0.0332
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

  [stress_xx_time]
    type = ElementAverageValue
    variable = stress_xx
  []

  [./computation_time]
    type = PerfGraphData
    section_name = Root
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
   [./fields]
    type = LineValueSampler
    num_points = 601
    start_point = '-4410 0 0'
    end_point = '4410 0 0'
    sort_by = id
    variable = 'eta x_Al u_x interface_energy'
   [../]

   [./derived_fields]
    type = LineValueSampler
    num_points = 601
    start_point = '-4410 0 0'
    end_point = '4410 0 0'
    sort_by = id
    variable = 'total_strain_xx stress_xx'
   [../]
[]

#[Adaptivity]
   # marker = 'errorfrac1'
   # max_h_level = 2
   # initial_steps = 5

   #[./Indicators]
   # [./error_eta]
   #     type = GradientJumpIndicator
   #     variable = eta
   # [../]
   #[../]

   #[./Markers]
   # [./errorfrac1]
   #     type = ErrorFractionMarker
   #     refine = 0.3
   #     coarsen = 0.1
   #     indicator = 'error_eta'
   # [../]
   #[../]
#[]


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
  l_tol = 1e-4
  nl_max_its = 15
  nl_rel_tol = 1.0e-8
  nl_abs_tol = 1.0e-10

  end_time = 4e9

  [./TimeStepper]
   type = SolutionTimeAdaptiveDT
   dt = 1.0
  [../]

  #[./TimeStepper]
  #  dt = 1.0
  #  type = IterationAdaptiveDT
  #  cutback_factor = 0.4
  #  growth_factor = 1.5
  #  optimal_iterations = 10
  #[../]
[]


[Outputs]

  file_base = q1D_twophase_NiAl_transient_ti
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
