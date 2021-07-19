#In this example, the mobility of liquid phase
#is 10e5 times that of Al in B2-NiAl phase
#This is a model Ni-Al alloy problem
#alpha phase- FCC_A1
#beta phase- FCC_L12
#gamma phase- BCC_B2
#delta phase- Al3NI2
#epsilon phase- Liquid
#T = 1000K
#Simulation details:
#Interfacial energy = 0.5 J/m2
#interfacial width = 1.5e-6 m
#Length = 300e-6 m
#characteristic_length = 0.05e-6 m
#charactersitic_time   =  324 s
#Nondimensional length = (Length/characteristic_length)
#No. of finite elements= 6.0(Length/interface_width)

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1200
  xmin = -3000
  xmax =  3000
  elem_type = EDGE2
[]


[Variables]

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

[./phi_delta]
  order = FIRST
  family = LAGRANGE
[../]

[./phi_epsilon]
  order= FIRST
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


[BCs]

  [./alpha_pf]
    type = NeumannBC
    variable = 'phi_alpha'
    boundary = 'left right'
    value = 0
  [../]
  
  [./beta_pf]
    type = NeumannBC
    variable = 'phi_beta'
    boundary = 'left right'
    value = 0
  [../]
  
   [./MassFluxBC]
    type = NeumannBC
    variable = 'mu_Al'
    boundary = 'left right'
    value = 0
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
  
  [./ic_phi_delta]
    type = FunctionIC
    variable = phi_delta
    function = func_phi_delta
  [../]
  
   [./ic_phi_epsilon]
    type = FunctionIC
    variable = phi_epsilon
    function = func_phi_epsilon
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
  
[]

[Functions]

  #Width of Ni3Al  = 20e-6 m
  #Width of NiAl   = 10e-6 m
  #Width of Al3Ni2 = 10e-6 m
  
  [./func_phi_alpha]
    type = ParsedFunction
    value = 'if(x<=-400, 1.0, 0)'
  [../]
  
  [./func_phi_beta]
    type = ParsedFunction
    value = 'if(x>-400 & x<=0.0, 1.0, 0)'
  [../]
  
  [./func_phi_gamma]
    type = ParsedFunction
    value = 'if(x>0 & x<=200, 1.0, 0)'
  [../]
  
  [./func_phi_delta]
    type = ParsedFunction
    value = 'if(x>200 & x<=400, 1.0, 0)'
  [../]
  
  [./func_phi_epsilon]
    type = ParsedFunction
    value = 'if(x>400, 1.0, 0)'
  [../]
  

  [./func_xAl]
    type = ParsedFunction
    vars = 'alpha beta gamma delta epsilon'
    vals = 'func_phi_alpha  func_phi_beta func_phi_gamma func_phi_delta func_phi_epsilon'
    value = '1.28E-1*alpha + 2.4E-1*beta + 5.41418E-1*gamma + 6.80224E-1*delta + 8.93E-1*epsilon'
  [../]
 
  [./func_muAl]
    type = ParsedFunction
    vars = 'alpha beta gamma delta epsilon'
    vals = 'func_phi_alpha  func_phi_beta func_phi_gamma func_phi_delta func_phi_epsilon'
    value = '-1.63E1*alpha -1.63E1*beta + 6.4E0*gamma + 1.625E1*delta + 1.6251E1*epsilon'
  [../]
  
[]



[Kernels]
  
  #Variable that this kernel operates on
  #is the diffusion potential of comp B
  
  [./PhaseConc]
    type = MultiPhaseConstraintMu
    variable = mu_Al
    xB = x_Al
    phase_alpha = phi_alpha
    phase_beta = phi_beta
    phase_gamma = phi_gamma
    phase_delta = phi_delta
    phase_epsilon = phi_epsilon
    xB_gamma = xB_gamma
    xB_delta = xB_delta
    xB_epsilon = xB_epsilon
    inv_B_tf_gamma = inv_B_tf_gamma
    inv_B_tf_delta = inv_B_tf_delta
    inv_B_tf_epsilon = inv_B_tf_epsilon
    h_gamma = h_gamma
    h_delta = h_delta
    h_epsilon = h_epsilon
  [../]

  ######################
  # Diffusion Equation #
  ######################
    
  [./Al_balance]
    type = BinaryMultiPhaseMassBalance
    variable = x_Al
    B_diff_pot = mu_Al
    phase_alpha = phi_alpha
    phase_beta = phi_beta
    phase_gamma = phi_gamma
    phase_delta = phi_delta
    phase_epsilon = phi_epsilon
    h_alpha = h_alpha
    h_beta = h_beta
    h_gamma = h_gamma
    h_delta = h_delta
    h_epsilon = h_epsilon
    inv_B_tf_gamma = inv_B_tf_gamma
    inv_B_tf_delta = inv_B_tf_delta
    inv_B_tf_epsilon = inv_B_tf_epsilon
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
    phase_4 = phi_delta
    phase_5 = phi_epsilon
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_alpha
    A_chem_pot_2 = A_chem_pot_beta
    xB_1 = xB_alpha
    xB_2 = xB_beta
    dh = dhbeta_dphialpha     #Required by the residual
    d2h = d2hbeta_dphialpha2  #Required by the diag jacobian
    d2h_2 = d2hbeta_dphialpha_dphibeta
    d2h_3 = d2hbeta_dphialpha_dphigamma
    d2h_4 = d2hbeta_dphialpha_dphidelta
    d2h_5 = d2hbeta_dphialpha_dphiepsilon
    mob_name = L_phi
    nd_factor = nd_factor
  [../]
  
  
  [./BinaryMultiPhaseDrivingForce_alphagamma]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_alpha     #Phase_1
    phase_2 = phi_gamma      #Coupled phase_2
    phase_3 = phi_beta       #Coupled phase_3
    phase_4 = phi_delta
    phase_5 = phi_epsilon
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_alpha
    A_chem_pot_2 = A_chem_pot_gamma
    xB_1 = xB_alpha
    xB_2 = xB_gamma
    dh = dhgamma_dphialpha     #Required by the residual
    d2h = d2hgamma_dphialpha2  #Required by the diag jacobian
    d2h_2 = d2hgamma_dphialpha_dphigamma
    d2h_3 = d2hgamma_dphialpha_dphibeta
    d2h_4 = d2hgamma_dphialpha_dphidelta
    d2h_5 = d2hgamma_dphialpha_dphiepsilon
    mob_name = L_phi
    nd_factor = nd_factor
  [../]
  
   [./BinaryMultiPhaseDrivingForce_alphadelta]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_alpha     #Phase_1
    phase_2 = phi_delta      #Coupled phase_2
    phase_3 = phi_beta       #Coupled phase_3
    phase_4 = phi_gamma
    phase_5 = phi_epsilon
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_alpha
    A_chem_pot_2 = A_chem_pot_delta
    xB_1 = xB_alpha
    xB_2 = xB_delta
    dh = dhdelta_dphialpha     #Required by the residual
    d2h = d2hdelta_dphialpha2  #Required by the diag jacobian
    d2h_2 = d2hdelta_dphialpha_dphidelta
    d2h_3 = d2hdelta_dphialpha_dphibeta
    d2h_4 = d2hdelta_dphialpha_dphigamma
    d2h_5 = d2hdelta_dphialpha_dphiepsilon
    mob_name = L_phi
    nd_factor = nd_factor
  [../]
  
  [./BinaryMultiPhaseDrivingForce_alphaepsilon]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_alpha     #Phase_1
    phase_2 = phi_epsilon    #Coupled phase_2
    phase_3 = phi_beta       #Coupled phase_3
    phase_4 = phi_gamma
    phase_5 = phi_delta
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_alpha
    A_chem_pot_2 = A_chem_pot_epsilon
    xB_1 = xB_alpha
    xB_2 = xB_epsilon
    dh = dhepsilon_dphialpha     #Required by the residual
    d2h = d2hepsilon_dphialpha2  #Required by the diag jacobian
    d2h_2 = d2hepsilon_dphialpha_dphiepsilon
    d2h_3 = d2hepsilon_dphialpha_dphibeta
    d2h_4 = d2hepsilon_dphialpha_dphigamma
    d2h_5 = d2hepsilon_dphialpha_dphidelta
    mob_name = L_phi
    nd_factor = nd_factor
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
    eta4 = phi_delta
    eta5 = phi_epsilon
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
    phase_4 = phi_delta
    phase_5 = phi_epsilon
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_beta
    A_chem_pot_2 = A_chem_pot_alpha
    xB_1 = xB_beta
    xB_2 = xB_alpha
    dh = dhalpha_dphibeta
    d2h = d2halpha_dphibeta2
    d2h_2 = d2halpha_dphibeta_dphialpha
    d2h_3 = d2halpha_dphibeta_dphigamma
    d2h_4 = d2halpha_dphibeta_dphidelta
    d2h_5 = d2halpha_dphibeta_dphiepsilon
    mob_name = L_phi
    nd_factor = nd_factor
  [../]
  
  [./BinaryMultiPhaseDrivingForce_betagamma]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_beta
    phase_2 = phi_gamma
    phase_3 = phi_alpha
    phase_4 = phi_delta
    phase_5 = phi_epsilon
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_beta
    A_chem_pot_2 = A_chem_pot_gamma
    xB_1 = xB_beta
    xB_2 = xB_gamma
    dh = dhgamma_dphibeta
    d2h = d2hgamma_dphibeta2
    d2h_2 = d2hgamma_dphibeta_dphigamma
    d2h_3 = d2hgamma_dphibeta_dphialpha
    d2h_4 = d2hgamma_dphibeta_dphidelta
    d2h_5 = d2hgamma_dphibeta_dphiepsilon
    mob_name = L_phi
    nd_factor = nd_factor
  [../]
  
   [./BinaryMultiPhaseDrivingForce_betadelta]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_beta
    phase_2 = phi_delta
    phase_3 = phi_gamma
    phase_4 = phi_alpha
    phase_5 = phi_epsilon
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_beta
    A_chem_pot_2 = A_chem_pot_delta
    xB_1 = xB_beta
    xB_2 = xB_delta
    dh = dhdelta_dphibeta
    d2h = d2hdelta_dphibeta2
    d2h_2 = d2hdelta_dphibeta_dphidelta
    d2h_3 = d2hdelta_dphibeta_dphigamma
    d2h_4 = d2hdelta_dphibeta_dphialpha
    d2h_5 = d2hdelta_dphibeta_dphiepsilon
    mob_name = L_phi
    nd_factor = nd_factor
  [../]
  
   [./BinaryMultiPhaseDrivingForce_betaepsilon]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_beta
    phase_2 = phi_epsilon
    phase_3 = phi_gamma
    phase_4 = phi_alpha
    phase_5 = phi_delta
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_beta
    A_chem_pot_2 = A_chem_pot_epsilon
    xB_1 = xB_beta
    xB_2 = xB_epsilon
    dh = dhepsilon_dphibeta
    d2h = d2hepsilon_dphibeta2
    d2h_2 = d2hepsilon_dphibeta_dphiepsilon
    d2h_3 = d2hepsilon_dphibeta_dphigamma
    d2h_4 = d2hepsilon_dphibeta_dphialpha
    d2h_5 = d2hepsilon_dphibeta_dphidelta
    mob_name = L_phi
    nd_factor = nd_factor
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
    eta4 = phi_delta
    eta5 = phi_epsilon
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
  
  [./BinaryMultiPhaseDrivingForce_gammaalpha]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_gamma
    phase_2 = phi_alpha
    phase_3 = phi_beta
    phase_4 = phi_delta
    phase_5 = phi_epsilon
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_gamma
    A_chem_pot_2 = A_chem_pot_alpha
    xB_1 = xB_gamma
    xB_2 = xB_alpha
    dh = dhalpha_dphigamma
    d2h = d2halpha_dphigamma2
    d2h_2 = d2halpha_dphigamma_dphialpha
    d2h_3 = d2halpha_dphigamma_dphibeta
    d2h_4 = d2halpha_dphigamma_dphidelta
    d2h_5 = d2halpha_dphigamma_dphiepsilon
    mob_name = L_phi
    nd_factor = nd_factor
  [../]
   
  [./BinaryMultiPhaseDrivingForce_gammabeta]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_gamma
    phase_2 = phi_beta
    phase_3 = phi_alpha
    phase_4 = phi_delta
    phase_5 = phi_epsilon
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_gamma
    A_chem_pot_2 = A_chem_pot_beta
    xB_1 = xB_gamma
    xB_2 = xB_beta
    dh = dhbeta_dphigamma
    d2h = d2hbeta_dphigamma2
    d2h_2 = d2hbeta_dphigamma_dphibeta
    d2h_3 = d2hbeta_dphigamma_dphialpha
    d2h_4 = d2hbeta_dphigamma_dphidelta
    d2h_5 = d2hbeta_dphigamma_dphiepsilon
    mob_name = L_phi
    nd_factor = nd_factor
  [../]
  
  [./BinaryMultiPhaseDrivingForce_gammadelta]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_gamma
    phase_2 = phi_delta
    phase_3 = phi_alpha
    phase_4 = phi_beta
    phase_5 =  phi_epsilon
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_gamma
    A_chem_pot_2 = A_chem_pot_delta
    xB_1 = xB_gamma
    xB_2 = xB_delta
    dh = dhdelta_dphigamma
    d2h = d2hdelta_dphigamma2
    d2h_2 = d2hdelta_dphigamma_dphidelta
    d2h_3 = d2hdelta_dphigamma_dphialpha
    d2h_4 = d2hdelta_dphigamma_dphibeta
    d2h_5 = d2hdelta_dphigamma_dphiepsilon
    mob_name = L_phi
    nd_factor = nd_factor
  [../]
  
  [./BinaryMultiPhaseDrivingForce_gammaepsilon]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_gamma
    phase_2 = phi_epsilon
    phase_3 = phi_alpha
    phase_4 = phi_beta
    phase_5 =  phi_delta
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_gamma
    A_chem_pot_2 = A_chem_pot_epsilon
    xB_1 = xB_gamma
    xB_2 = xB_epsilon
    dh = dhepsilon_dphigamma
    d2h = d2hepsilon_dphigamma2
    d2h_2 = d2hepsilon_dphigamma_dphiepsilon
    d2h_3 = d2hepsilon_dphigamma_dphialpha
    d2h_4 = d2hepsilon_dphigamma_dphibeta
    d2h_5 = d2hepsilon_dphigamma_dphidelta
    mob_name = L_phi
    nd_factor = nd_factor
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
    eta4 = phi_delta
    eta5 = phi_epsilon
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
  
  #######################
  # Allen-Cahn equation #
  # for phase-delta     #
  #######################
   
  [./BinaryMultiPhaseDrivingForce_deltaalpha]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_delta
    phase_2 = phi_alpha
    phase_3 = phi_beta
    phase_4 = phi_gamma
    phase_5 = phi_epsilon
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_delta
    A_chem_pot_2 = A_chem_pot_alpha
    xB_1 = xB_delta
    xB_2 = xB_alpha
    dh = dhalpha_dphidelta
    d2h = d2halpha_dphidelta2
    d2h_2 = d2halpha_dphidelta_dphialpha
    d2h_3 = d2halpha_dphidelta_dphibeta
    d2h_4 = d2halpha_dphidelta_dphigamma
    d2h_5 = d2halpha_dphidelta_dphiepsilon
    mob_name = L_phi
    nd_factor = nd_factor
  [../]
  
  [./BinaryMultiPhaseDrivingForce_deltabeta]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_delta
    phase_2 = phi_beta
    phase_3 = phi_gamma
    phase_4 = phi_alpha
    phase_5 = phi_epsilon
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_delta
    A_chem_pot_2 = A_chem_pot_beta
    xB_1 = xB_delta
    xB_2 = xB_beta
    dh = dhbeta_dphidelta
    d2h = d2hbeta_dphidelta2
    d2h_2 = d2hbeta_dphidelta_dphibeta
    d2h_3 = d2hbeta_dphidelta_dphigamma
    d2h_4 = d2hbeta_dphidelta_dphialpha
    d2h_5 = d2hbeta_dphidelta_dphiepsilon
    mob_name = L_phi
    nd_factor = nd_factor
  [../]
  
  [./BinaryMultiPhaseDrivingForce_deltagamma]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_delta
    phase_2 = phi_gamma
    phase_3 = phi_alpha
    phase_4 = phi_beta
    phase_5 = phi_epsilon
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_delta
    A_chem_pot_2 = A_chem_pot_gamma
    xB_1 = xB_delta
    xB_2 = xB_gamma
    dh = dhgamma_dphidelta
    d2h = d2hgamma_dphidelta2
    d2h_2 = d2hgamma_dphidelta_dphigamma
    d2h_3 = d2hgamma_dphidelta_dphialpha
    d2h_4 = d2hgamma_dphidelta_dphibeta
    d2h_5 = d2hgamma_dphidelta_dphiepsilon
    mob_name = L_phi
    nd_factor = nd_factor
  [../]
  
   [./BinaryMultiPhaseDrivingForce_deltaepsilon]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_delta
    phase_2 = phi_epsilon
    phase_3 = phi_alpha
    phase_4 = phi_beta
    phase_5 = phi_gamma
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_delta
    A_chem_pot_2 = A_chem_pot_epsilon
    xB_1 = xB_delta
    xB_2 = xB_epsilon
    dh = dhepsilon_dphidelta
    d2h = d2hepsilon_dphidelta2
    d2h_2 = d2hepsilon_dphidelta_dphiepsilon
    d2h_3 = d2hepsilon_dphidelta_dphialpha
    d2h_4 = d2hepsilon_dphidelta_dphibeta
    d2h_5 = d2hepsilon_dphidelta_dphigamma
    mob_name = L_phi
    nd_factor = nd_factor
  [../]
  
  #This kernel implements the double well
  #The form of the double-well barrier is 
  #f(phi) = phi^(2)(1 - phi^(2))
  #The free energy is explicitly coded
  #within the kernel
  
   [./doublewell_phi_delta]
    type = MultiPhaseDoubleWell
    variable = phi_delta
    eta2 = phi_alpha
    eta3 = phi_beta
    eta4 = phi_gamma 
    eta5 = phi_epsilon
    gamma = 1.5
    mob_name = L_phi
    barrier_height = m
  [../]

  [./Curvature_delta]
    type = ACInterface
    variable = phi_delta
    kappa_name = kappa
    mob_name = L_phi
  [../]
  
  [./dphi_delta_dt]
    type = TimeDerivative
    variable = phi_delta
  [../]
  
    
  #######################
  # Allen-Cahn equation #
  # for phase-epsilon   #
  #######################
   
  [./BinaryMultiPhaseDrivingForce_epsilonalpha]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_epsilon
    phase_2 = phi_alpha
    phase_3 = phi_beta
    phase_4 = phi_gamma
    phase_5 = phi_delta
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_epsilon
    A_chem_pot_2 = A_chem_pot_alpha
    xB_1 = xB_epsilon
    xB_2 = xB_alpha
    dh = dhalpha_dphiepsilon
    d2h = d2halpha_dphiepsilon2
    d2h_2 = d2halpha_dphiepsilon_dphialpha
    d2h_3 = d2halpha_dphiepsilon_dphibeta
    d2h_4 = d2halpha_dphiepsilon_dphigamma
    d2h_5 = d2halpha_dphiepsilon_dphidelta
    mob_name = L_phi
    nd_factor = nd_factor
  [../]
  
  [./BinaryMultiPhaseDrivingForce_epsilonbeta]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_epsilon
    phase_2 = phi_beta
    phase_3 = phi_gamma
    phase_4 = phi_alpha
    phase_5 = phi_delta
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_epsilon
    A_chem_pot_2 = A_chem_pot_beta
    xB_1 = xB_epsilon
    xB_2 = xB_beta
    dh = dhbeta_dphiepsilon
    d2h = d2hbeta_dphiepsilon2
    d2h_2 = d2hbeta_dphiepsilon_dphibeta
    d2h_3 = d2hbeta_dphiepsilon_dphigamma
    d2h_4 = d2hbeta_dphiepsilon_dphialpha
    d2h_5 = d2hbeta_dphiepsilon_dphidelta
    mob_name = L_phi
    nd_factor = nd_factor
  [../]
  
  [./BinaryMultiPhaseDrivingForce_epsilongamma]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_epsilon
    phase_2 = phi_gamma
    phase_3 = phi_alpha
    phase_4 = phi_beta
    phase_5 = phi_delta
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_epsilon
    A_chem_pot_2 = A_chem_pot_gamma
    xB_1 = xB_epsilon
    xB_2 = xB_gamma
    dh = dhgamma_dphiepsilon
    d2h = d2hgamma_dphiepsilon2
    d2h_2 = d2hgamma_dphiepsilon_dphigamma
    d2h_3 = d2hgamma_dphiepsilon_dphialpha
    d2h_4 = d2hgamma_dphiepsilon_dphibeta
    d2h_5 = d2hgamma_dphiepsilon_dphidelta
    mob_name = L_phi
    nd_factor = nd_factor
  [../]
  
   [./BinaryMultiPhaseDrivingForce_epsilondelta]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_epsilon
    phase_2 = phi_delta
    phase_3 = phi_alpha
    phase_4 = phi_beta
    phase_5 = phi_gamma
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_epsilon
    A_chem_pot_2 = A_chem_pot_delta
    xB_1 = xB_epsilon
    xB_2 = xB_delta
    dh = dhdelta_dphiepsilon
    d2h = d2hdelta_dphiepsilon2
    d2h_2 = d2hdelta_dphiepsilon_dphidelta
    d2h_3 = d2hdelta_dphiepsilon_dphialpha
    d2h_4 = d2hdelta_dphiepsilon_dphibeta
    d2h_5 = d2hdelta_dphiepsilon_dphigamma
    mob_name = L_phi
    nd_factor = nd_factor
  [../]
  
  #This kernel implements the double well
  #The form of the double-well barrier is 
  #f(phi) = phi^(2)(1 - phi^(2))
  #The free energy is explicitly coded
  #within the kernel
  
   [./doublewell_phi_epsilon]
    type = MultiPhaseDoubleWell
    variable = phi_epsilon
    eta2 = phi_alpha
    eta3 = phi_beta
    eta4 = phi_gamma 
    eta5 = phi_delta
    gamma = 1.5
    mob_name = L_phi
    barrier_height = m
  [../]

  [./Curvature_epsilon]
    type = ACInterface
    variable = phi_epsilon
    kappa_name = kappa
    mob_name = L_phi
  [../]
  
  [./dphi_epsilon_dt]
    type = TimeDerivative
    variable = phi_epsilon
  [../]
  
[]

[Materials]
  ##Phenomenological coeff###
  
  #nd_factor = R*T/(Vm*m)
  #The dimensional value of LB(B2,Al,Al,Ni) = 1.2515e-23 molm2/Js
  #characteristic_Onsager_constant Lc = 9.2808e-22 molm2/Js
  
    [./ConstantFieldProperties]
     type = GenericConstantMaterial
     prop_names =  'L_phi       kappa   m   nd_factor L_BB_gamma   L_BB_delta  L_BB_epsilon'
     prop_values = '2.2772e-4   112.5  1.0  55.4267   1.35e-2        1.35e-2     1.35e3'
  [../]
  
   [./ZeroFieldProperties]
     type = GenericConstantMaterial
     prop_names =  'dL_BB_muB_gamma dL_BB_muB_delta dL_BB_muB_epsilon' 
     prop_values = '0 0 0'
  [../] 
  
  [./InterpolationFunction]
    type = QuantInterpolationFunction
    phase_alpha = phi_alpha
    phase_beta = phi_beta
    phase_gamma = phi_gamma
    phase_delta = phi_delta
    phase_epsilon = phi_epsilon
  [../]  
  
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
  
   [./BinaryConjugatePhaseMaterial_Al3Ni2]
    type = BinaryConjugatePhaseMaterial
    A_chem_pot = A_chem_pot_delta
    xB = xB_delta
    inv_B_tf = inv_B_tf_delta
    B_diff_pot = mu_Al
    molar_volume = 1.0
    table_object = table_Al3Ni2_data
    #outputs = exodus
  [../]
  
  [./BinaryConjugatePhaseMaterial_Liq]
    type = BinaryConjugatePhaseMaterial
    A_chem_pot = A_chem_pot_epsilon
    xB = xB_epsilon
    inv_B_tf = inv_B_tf_epsilon
    B_diff_pot = mu_Al
    molar_volume = 1.0
    table_object = table_Liquid_data
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
  
  #[./BinaryConjugateKineticMaterial_B2]
    #type = BinaryConjugateKineticMaterial
    #L_BB = L_BB_gamma_tc
    #dL_BB_muB = dL_BB_muB_gamma_tc
    #B_diff_pot = mu_Al
    #table_object = table_B2_mobility_data
    #outputs = exodus
  #[../]
  
   #[./BinaryConjugateKineticMaterial_Liquid]
   # type = BinaryConjugateKineticMaterial
   # L_BB = L_BB_epsilon_tc
   # dL_BB_muB = dL_BB_muB_epsilon_tc
   # B_diff_pot = mu_Al
   # table_object = table_Liquid_mobility_data
   # outputs = exodus
   #[../]
   
[]

[UserObjects]

  #This userobjcts reads the data from the table 
  #for a given diffusion potential
  #and provides the interpolated grand-potential,
  #mole fraction and inverse of the thermodynamic_factor
  #to the material TabulatedPhaseMaterial
  
  [./table_FCC_data]
    type = BinaryConjugatePhaseData
    table_name = NiAl_data/NiAl_invert_data_FCC.csv
  [../] 
  
  [./table_L12_data]
    type = BinaryConjugatePhaseData
    table_name = NiAl_data/NiAl_invert_data_L12.csv
  [../]
  
   [./table_B2_data]
    type = BinaryConjugatePhaseData
    table_name = NiAl_data/NiAl_invert_data_B2.csv
  [../] 
  
  [./table_Al3Ni2_data]
    type = BinaryConjugatePhaseData
    table_name = NiAl_data/NiAl_invert_data_Al3Ni2.csv
  [../]
  
  [./table_Liquid_data]
    type = BinaryConjugatePhaseData
    table_name = NiAl_data/NiAl_invert_data_Liquid.csv
  [../]
  
  #Kinetic data files
  
  [./table_FCC_mobility_data]
    type = BinaryConjugateMobilityData
    table_name = NiAl_data/NiAl_mobility_data_FCC_invert.csv
  [../]
  
  [./table_L12_mobility_data]
   type = BinaryConjugateMobilityData
   table_name = NiAl_data/NiAl_mobility_data_L12_invert.csv
  [../]
  
  #[./table_B2_mobility_data]
  #  type = BinaryConjugateMobilityData
  #  table_name = NiAl_data/NiAl_mobility_data_B2_invert.csv
  #[../]
  
  #[./table_Liquid_mobility_data]
  # type = BinaryConjugateMobilityData
  # table_name = NiAl_data/NiAl_mobility_data_Liquid_invert.csv
  #[../]
  
  #[./table_Al3Ni2_mobility_data]
  #  type = BinaryConjugateMobilityData
  #  table_name = NiAl_data/NiAl_mobility_data_Al3Ni2_invert.csv
  #[../]
  
[]

[Postprocessors]
 
  [./computation_time]
    type = PerfGraphData
    section_name = Transient::PicardSolve
    data_type = total
    execute_on = timestep_end
  [../]

  [./pos_A1_L12]
    type = FindValueOnLine
    target = 0.5
    v = phi_alpha
    start_point = '-3000 0 0'
    end_point = '3000 0 0'
    tol = 1e-8
    execute_on = 'initial timestep_end'
  [../]

  [./pos_L12_B2]
    type = FindValueOnLine
    target = .5
    v = phi_gamma
    start_point = '-3000 0 0'
    end_point = '20 0 0'
    tol = 1e-8
    execute_on = 'initial timestep_end'
  [../]
 
  [./pos_B2_Al3Ni2]
   type = FindValueOnLine
   target = .5
   v = phi_gamma
   start_point = '20 0 0'
   end_point = '3000 0 0'
   tol = 1e-8
   execute_on = 'initial timestep_end'
 [../]
 
 [./pos_Al3Ni2_liq]
  type = FindValueOnLine
  target = .5
  v = phi_epsilon
  start_point = '-3000 0 0'
  end_point = '3000 0 0'
  tol = 1e-8
  execute_on = 'initial timestep_end'
 [../] 

  [./volumefraction_alpha]
    type = VolumeFraction
    variable = phi_alpha
    execute_on = 'initial timestep_end'
  [../]
  
  [./volumefraction_beta]
    type = VolumeFraction
    variable = phi_beta
    execute_on = 'initial timestep_end'
 [../]
 
 [./volumefraction_gamma]
  type = VolumeFraction
  variable = phi_gamma
  execute_on = 'initial timestep_end'
 [../]

  [./volumefraction_delta]
    type = VolumeFraction
    variable = phi_delta
    execute_on = 'initial timestep_end'
  [../]
  
  [./volumefraction_epsilon]
    type = VolumeFraction
    variable = phi_epsilon
    execute_on = 'initial timestep_end'
 [../]

[]

#[VectorPostprocessors]

  #[./xAl]
  #  type = LineValueSampler
  #  start_point = '-600 0 0'
  #  end_point = '600 0 0'
  #  num_points = 1200
  #  variable = 'x_Al'
  #  sort_by = id
  #  execute_on = 'initial timestep_end'
  #[../]


#[]

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
  petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm       ilu          nonzero'
  
  l_tol = 1.0e-5
  l_max_its = 15
  nl_rel_tol = 1.0e-8
  nl_abs_tol = 1.0e-10
  nl_max_its = 30
  
  #time = 4.444e2
  #real-time = 1000hr
  
  #1000 hr = 1.11111111e4
   
  end_time = 1.7e6 #non-dimensional time = 10K hrs
  
  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    percent_change = 0.1
    dt = 1.0  #real_dt = 324sec
  [../]
  
[]

[Outputs]

  [./csv]
    type = CSV
    execute_on = 'initial timestep_end'
  [../]
  
   [pgraph]
    type = PerfGraphOutput
    execute_on = 'final'
    level = 1  
    heaviest_branch = false
  []

  file_base = NiAl_diffusion_couple_ti_fast
  exodus = true
  
[]

[Debug]
  show_var_residual_norms = true
[] 
