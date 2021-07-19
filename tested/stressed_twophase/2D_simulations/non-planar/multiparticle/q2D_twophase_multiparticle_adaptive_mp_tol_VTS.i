#This test problem combines VTS
#with grand-potential model
#for a multi-particle problem
#T = 1473 K
#Al_diff_pot = -1.35156e5 J/mol
#interface_width = 50e-9 m
#interface_energy = 36.2e-3 J/m2
#for cubic materials
#The characteristic length and time
# are 1.1338E-9 m and 4.8411E-6 s

#Dimensional values:
#Length = 5e-6 m
#Mean radius = 1.0e-7 m

#Non-dimensional values:
#Length = 4410
#Mean radius = 88.2
#Number of partcles = 160


[Mesh]
  type = GeneratedMesh
  dim  = 2
   nx  = 150
   ny  = 150
 xmin  = -2205
 xmax  =  2205
 ymin  = -2205
 ymax  =  2205
 elem_type = QUAD4
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

  [./u_y]
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

  [./ic_xAl]
   type = FunctionIC
   variable = x_Al
   function = ic_xAl_fun
  [../]

  #[./ic_eta]
  #  type = SmoothCircleIC
  #  variable = eta
  #  x1 = 0
  #  y1 = 0
  #  radius = 176.4
  #  invalue = 1
  #  outvalue = 0
  #  int_width = 0
  #[../]

  #[./ic_xAl]
  # type = SmoothCircleIC
  # variable = x_Al
  # x1 = 0
  # y1 = 0
  # radius = 176.4
  # invalue = 2.34E-1
  # outvalue = 1.9E-1
  # int_width = 0
  # #function = ic_func_xAl
  #[../]

[]

[Functions]
  [./ic_eta_fun]
    type = CircleFunctionIC
    inside_value = 1.0
    outside_value = 0.0
    table_name = 'outputfile_mp.dat'
  [../]

  [./ic_xAl_fun]
    type = ParsedFunction
    vars = 'eta  xAl_beta xAl_alpha'
    vals = 'ic_eta_fun 2.34E-1 1.92E-1'
    value = 'xAl_beta*eta + xAl_alpha*(1-eta)'
  [../]
[]

[AuxVariables]

  [./delta_omega]
   order = CONSTANT
   family = MONOMIAL
  [../]

  [./grad_eta_0]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./strain_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./strain_xy]
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

  [./ChemicalPotential]
   type = ChemicalPotential
   variable = delta_omega
   dh = dh
  [../]

  [./stress_xx]
    type = RankTwoAux
    variable = stress_xx
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
  [../]

  [./stress_yy]
    type = RankTwoAux
    variable = stress_yy
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
  [../]

  [./stress_xy]
    type = RankTwoAux
    variable = stress_xy
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
  [../]

  #alpha_total_strain = beta_total_strain

  [./strain_xx]
    type = RankTwoAux
    variable = strain_xx
    rank_two_tensor = alpha_total_strain
    index_i = 0
    index_j = 0
  [../]

  [./strain_yy]
    type = RankTwoAux
    variable = strain_yy
    rank_two_tensor = alpha_total_strain
    index_i = 1
    index_j = 1
  [../]

   [./strain_xy]
    type = RankTwoAux
    variable = strain_xy
    rank_two_tensor = alpha_total_strain
    index_i = 0
    index_j = 1
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

 #Due to symmetric BC
 #left_disp_x = 0
 #bottom_disp_y = 0

 #[./left_disp_x]
 # type = DirichletBC
 # variable = 'u_x'
 # boundary = 'left'
 # value = 0.0
 #[../]

 #[./bottom_disp_y]
 #  type = DirichletBC
 #  variable = 'u_y'
 #   boundary = 'bottom'
 #  value = 0
 #[../]


 [./Periodic]
   [./auto]
     variable = 'eta mu_Al x_Al u_x u_y'
     auto_direction = 'x y'
   [../]
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
  #  nd_factor = ndf_el
  #[../]
[]


[Materials]

  [./consts]
    type = GenericConstantMaterial
    prop_names  = 'L_phi     kappa     m    nd_factor  ndf_el'
    prop_values = '0.0069    112.5     1.0   12.7854  6.5804E3'
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

  #[./beta_prefactor]
  #  type =  EigenStrainPrefactor
  #  base_name = beta
  #  R    = 8.314
  #  T    = 1473
  #  overall_xB = 0.19
  #  xB_eqm     = 2.3073E-1
  #  B_tf_eqm   = 2.86035E5
  #  B_diff_pot_eqm = -1.08531E5
  #  B_diff_pot     =  mu_Al
  #[../]

  [./beta_eigen]
    type = ComputeEigenstrain
    base_name = beta
    eigenstrain_name = eigen
    #prefactor = beta_prefactor
    eigen_base = '-0.003 0.0   0.0
                   0.0  -0.003 0.0
                   0.0    0.0  0.0'
    #outputs = exodus
  [../]

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

  [./TaylorVoigtHomogenization]
    type = TaylorVoigtHomogenization
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

  #[./ChemicalProperties_L12]
  # type = StrainDependentTaylorApproximation
  # base_name = beta
  # char_energy = 1.225e4
  # xB_eqm      = 2.3073E-1
  # B_tf_eqm    = 2.86035E5
  # B_diff_pot_eqm = -1.08531E5
  # A_chem_pot_eqm = -8.56263E4
  # B_diff_pot   = mu_Al
  #[../]

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
    L_BB_val = 0.0312
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

  [./VolumeFraction]
    type = VolumeFraction
    variable = eta
  [../]

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

  #[stress_xx_time]
  #  type = ElementAverageValue
  #  variable = stress_xx
  #[]

  #[strain_xx_avg]
  #  type = ElementAverageValue
  #  variable = total_strain_xx
  #[]

  [./computation_time]
    type = PerfGraphData
    section_name = Transient::PicardSolve
    data_type = total
    execute_on = 'timestep_end'
  [../]

  #[./position]
  #  type = FindValueOnLine
  #  target = 0.5
  #  v = eta
  #  start_point = '0      0 0'
  #  end_point   = '4410   0 0'
  #  tol = 1e-8
  #  execute_on = 'timestep_end'
  #[../]

  [./flood_count]
    type = FeatureFloodCount
    variable = eta
    compute_var_to_feature_map = true
    threshold = 0.5
    use_less_than_threshold_comparison = true
    execute_on = 'initial timestep_end'
  [../]

  [./mean_grain_area]
    type = AverageGrainVolume
    feature_counter = flood_count
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
    type = NodalValueSampler
    sort_by = id
    variable = 'eta x_Al u_x u_y'
    execute_on = 'initial timestep_end final'
   [../]

   [./derived_fields]
    type = LineValueSampler
    num_points = 601
    start_point = '0 0 0'
    end_point = '2205 0 0'
    sort_by = id
    variable = 'eta x_Al u_x u_y strain_xx strain_yy strain_xy stress_xx stress_yy stress_xy'
   [../]
[]

[Adaptivity]
    marker = 'errorfrac1'
    max_h_level = 2
    initial_steps = 4

   [./Indicators]
     [./error_xAl]
        type = GradientJumpIndicator
        variable = x_Al
    [../]
   [../]

   [./Markers]
    [./errorfrac1]
        type = ErrorFractionMarker
        refine = 0.3
        coarsen = 0.1
        indicator = 'error_xAl'
    [../]
   [../]
[]


[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]


[Executioner]
  type = Transient
  #scheme = bdf2
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type  -sub_pc_type '
  petsc_options_value = 'asm       lu'

  #petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type'
  #petsc_options_value = 'asm      ilu          nonzero'

  petsc_options = '-snes_converged_reason -ksp_converged_reason'
  #petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_mat_solver_package'
  #petsc_options_value = 'lu NONZERO superlu_dist'

  #petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  #petsc_options_value = 'lu strumpack'

  l_max_its = 30
  dtmin = 1e-8
  l_tol = 1e-6
  nl_max_its = 15
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8

  #1sec = 2.06564624E5

  start_time = 0
  #end_time = 2.06564624E7
  end_time = 4e9

  #num_steps = 30

  #[./TimeStepper]
  #type = SolutionTimeAdaptiveDT
  # dt = 1.0
  #[../]

  [./TimeStepper]
    dt = 1.0
    type = IterationAdaptiveDT
    cutback_factor = 0.2
    growth_factor = 1.3
    optimal_iterations = 12
    iteration_window = 2
  [../]
[]


[Outputs]

  file_base = q2D_twophase_multiparticle_adaptive_mp_tol_VTS
  #exodus = true

  [./exodus]
    type = Exodus
    execute_on = 'initial timestep_end final'
    interval = 50
  [../]

  [./out]
   type = Checkpoint
   interval = 50
   num_files= 2
  [../]

  [./csv]
    type = CSV
    execute_on = 'initial final'
    execute_postprocessors_on = ' initial timestep_end final'
    interval = 1
  [../]

  [pgraph]
   type = PerfGraphOutput
   execute_on = 'final'
   level = 1
   heaviest_branch = false
  []

[]
