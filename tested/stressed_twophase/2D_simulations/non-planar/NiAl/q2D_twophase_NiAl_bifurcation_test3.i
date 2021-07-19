#This test problem combines PRH
#with grand-potential model
#grand-potential model
#It reduces the domain using symmetry bc
#T = 1473 K
#Outer_radius = 5e-6 m
#Al_diff_pot = -1.35156e5 J/mol
#interface_width = 25e-9 m
#interface_energy = 36.2e-3 J/m2
#For a non-planar problem
#Radius of the precipitate is 0.1e-6 m (1/50)th of the outer radius
#Note that radii_start*num_circles = 5e-6m
#Nondimensional outer radius = 4410(= 5e-6 m)
#Non-dimensional ppt radius = 88.2
#Initial volume fraction = 0.04 %

[Mesh]
  [concentric_ring]
    type = MultiRadiiMeshGenerator
    num_sectors = 10
    radii = '3.675'         #I assume that this is the radius
                            #of the innermost circle
    rings = '1 '
    num_circles = 1200
    has_outer_square = false
    preserve_volumes = false
    portion = top_right
  []
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
    type = SmoothCircleIC
    variable = eta
    x1 = 0
    y1 = 0
    radius = 88.2
    invalue = 1
    outvalue = 0
    int_width = 0
  [../]

  [./ic_xAl]
   type = SmoothCircleIC
   variable = x_Al
   x1 = 0
   y1 = 0
   radius = 88.2
   invalue = 2.34E-1
   outvalue = 1.9E-1
   int_width = 0
   #function = ic_func_xAl
  [../]

[]

#[Functions]
#  [./ic_eta_fun]
#    type = ParsedFunction
#    vars = 'x1 y1 r'
#    vals = '0 0  1323'
#    value = 'if((x-x1)^2 + (y-y1)^2 <r*r,1,0)' 
#  [../]
  
#  [./ic_func_xAl]
#    type = ParsedFunction
#    vars = 'x1 y1 r   xAl_beta xAl_alpha'
#    vals = '0 0  1323  2.324E-1 1.9E-1'
#    value = 'if( (x-x1)^2 + (y-y1)^2 <= r*r, xAl_beta, xAl_alpha)'
#  [../]
#[]

[AuxVariables]

  [./delta_omega]
   order = CONSTANT
   family = MONOMIAL
  [../]

  [./grad_eta_0]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./radial_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./hoop_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./shear_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./radial_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./hoop_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]
    
  [./shear_strain]
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

  [./radial_stress]
    type = CylindricalRankTwoAux
    variable = radial_stress
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
    center_point = '0 0 0'
  [../]

  [./hoop_stress]
    type = CylindricalRankTwoAux
    variable = hoop_stress
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    center_point = '0 0 0'
  [../]

  [./shear_stress]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    variable = shear_stress
  [../]

  #alpha_total_strain = beta_total_strain

  [./radial_strain]
    type = CylindricalRankTwoAux
    variable = radial_strain
    rank_two_tensor = alpha_total_strain
    index_i = 0
    index_j = 0
    center_point = '0 0 0'
  [../]

  [./hoop_strain]
    type = CylindricalRankTwoAux
    variable = hoop_strain
    rank_two_tensor = alpha_total_strain
    index_i = 1
    index_j = 1
    center_point = '0 0 0'
  [../]

   [./shear_strain]
    type = RankTwoAux
    rank_two_tensor = alpha_total_strain
    index_i = 0
    index_j = 1
    variable = shear_strain
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

 [./left_disp_x]
  type = DirichletBC
  variable = 'u_x'
  boundary = 'left'
  value = 0.0
 [../] 

 [./bottom_disp_y]
   type = DirichletBC
   variable = 'u_y'
   boundary = 'bottom'
   value = 0
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
  
  [./QDrivingForceGradPhi3D]
    type = QDrivingForceGradPhi3D
    variable = eta
    mob_name = L_phi
    nd_factor = ndf_el
  [../]
[]


[Materials]

  [./consts]
    type = GenericConstantMaterial
    prop_names  = 'L_phi      kappa      m    nd_factor  ndf_el'
    prop_values = '0.0553    28.1250    1.0   6.3927    2.1683e3'
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

  [./drivingforce_gradphi]
    type = DrivingForceGradPhi
    eta  = eta
    strain_jump_name = strain_jump
  [../]
  
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

  #Assuming inhomogenous elasticity
  #The young's modulus of gamma_prime 
  #is equal to gamma
  
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
    eigen_base = '-0.003  0.0    0.0
                   0.0   -0.003  0.0
                   0.0   0.0     0.0'
    #outputs = exodus
  [../]
 
  ## The precipitate phase beta is stiffer than the matrix
  #For this reason E_beta > E_alpha
  
   [./beta_elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    base_name = beta
    youngs_modulus = 2.8528
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
  
  [./position]
    type = FindValueOnLine
    target = 0.5
    v = eta
    start_point = '0         0       0'
    end_point   = '3.1183e3 3.1183e3 0'
    tol = 1e-8
    execute_on = 'timestep_end'
  [../]

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
   [../]

   [./derived_fields]
    type = LineValueSampler
    num_points = 1201
    start_point = '0 0 0'
    end_point = '3.1183e3 3.1183e3 0'
    sort_by = id
    variable = 'eta x_Al u_x u_y radial_strain hoop_strain shear_strain radial_stress hoop_stress shear_stress'
   [../]
[]

#[Adaptivity]
   # marker = 'errorfrac1'
   # max_h_level = 2
   # initial_steps = 1
    
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
  #scheme = bdf2
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
  dtmin = 1e-8
  l_tol = 1e-4
  nl_max_its = 15
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  #1sec = 2.06564624E5

  start_time = 0
  end_time = 2.06564624E7
  #end_time = 4e5

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

  file_base = q2D_twophase_NiAl_bifurcation_test3
  #exodus = true

  [./exodus]
    type = Exodus
    execute_on = 'initial timestep_end final'
    interval = 50
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
