# ------------------------------------
# Secure the array for all cells
# ------------------------------------
function common_allocation(cellxmax, cellymax, nval)
    # define at cell center
    Qbase  = zeros(cellxmax, cellymax, nval)           # primitive variables
    volume = zeros(cellxmax, cellymax)                 # volume

    cellcenter = zeros(cellxmax, cellymax, 2)          # cell center
    
    wally = zeros(cellxmax, cellymax)                  # dis from wall

    # define at cell boundaries
    dx = zeros(cellxmax+1, cellymax)                   # distance from cell center
    dy = zeros(cellxmax, cellymax+1)                   # distance from cell center
    
    return Qbase, volume, cellcenter, wally, dx, dy
end 

# ------------------------------------
# Secure the array for mpi
# ------------------------------------
function common_mpi_allocation(cellxmax, cellymax, nval)
    # define at cell center
    Qbase  = zeros(cellxmax, cellymax, nval)           # primitive variables
    #Qbase_ave = zeros(cellxmax, cellymax, nval)        # average primitive variables 
    volume = zeros(cellxmax, cellymax)                 # volume

    cellcenter = zeros(cellxmax, cellymax, 2)          # cell center
    
    wally = zeros(cellxmax, cellymax)                  # dis from wall
    yplus = ones(cellxmax, cellymax) * 1e5             # yplus
    
    Qcon     = zeros(cellxmax, cellymax, nval)         # conserved variables
    Qcon_hat = zeros(cellxmax, cellymax, nval)         # conserved variables in the general coordinate system

    mu     = zeros(cellxmax, cellymax)                 # viscosity
    lambda = zeros(cellxmax, cellymax)                 # thermal Conductivity

    RHS = zeros(cellxmax, cellymax, nval)              # right hand side
    
    # define at cell boundaries
    dx = zeros(cellxmax+1, cellymax)                   # distance from cell center
    dy = zeros(cellxmax, cellymax+1)                   # distance from cell center
    
    vecAx = zeros(cellxmax+1, cellymax)                # vecter Ax
    vecAy = zeros(cellxmax, cellymax+1)                # vecter Ay

    E_adv_hat = zeros(cellxmax+1,   cellymax, nval)    # flux of advection in the x-direction
    F_adv_hat = zeros(  cellxmax, cellymax+1, nval)    # flux of advection in the y-direction

    E_vis_hat = zeros(cellxmax+1,   cellymax, nval)    # flux of viscosity in the x-direction
    F_vis_hat = zeros(  cellxmax, cellymax+1, nval)    # flux of viscosity in the y-direction

    return Qbase, volume, cellcenter, wally, yplus, dx, dy, Qcon, Qcon_hat, mu, lambda, 
            vecAx, vecAy, E_adv_hat, F_adv_hat, E_vis_hat, F_vis_hat, RHS
end 

function allocation_implicit(cellxmax, cellymax, nval)
    # define at cell center
    Qbasen     = zeros(cellxmax, cellymax, nval)        # primitive variables for inner iteration
    Qconn      = zeros(cellxmax, cellymax, nval)        # conserved variables for inner iteration
    Qconn_hat  = zeros(cellxmax, cellymax, nval)        # conserved variables in general coordinate system for inner iteration
    Qbasem     = zeros(cellxmax, cellymax, nval)        # primitive variables for inner iteration

    dtau         = zeros(cellxmax, cellymax)            # computational time steps

    A_adv_hat_p = zeros(cellxmax, cellymax, nval, nval) # Jacobian matrix A+ for one-wave approximation
    A_adv_hat_m = zeros(cellxmax, cellymax, nval, nval) # Jacobian matrix A- for one-wave approximation
    B_adv_hat_p = zeros(cellxmax, cellymax, nval, nval) # Jacobian matrix B+ for one-wave approximation
    B_adv_hat_m = zeros(cellxmax, cellymax, nval, nval) # Jacobian matrix B- for one-wave approximation
    A_beta_shig = zeros(cellxmax, cellymax)             # A sigma for one-wave approximation
    B_beta_shig = zeros(cellxmax, cellymax)             # B sigma for one-wave approximation

    jalphaP = zeros(cellxmax, cellymax)                 # Jacobian matrix for viscosity
    jbetaP  = zeros(cellxmax, cellymax)                 # Jacobian matrix for viscosity

    delta_Q      = zeros(cellxmax, cellymax, nval)      # delta Q for lusgs
    delta_Q_temp = zeros(cellxmax, cellymax, nval)      # temporary delta Q for lusgs

    D  = zeros(cellxmax, cellymax)                      # Diagonal of LHS
    Lx = zeros(cellxmax, cellymax, nval, nval)          # lower of LHS
    Ly = zeros(cellxmax, cellymax, nval, nval)          # lower of LHS
    Ux = zeros(cellxmax, cellymax, nval, nval)          # upper of LHS
    Uy = zeros(cellxmax, cellymax, nval, nval)          # upper of LHS

    LdQ = zeros(cellxmax, cellymax, nval)               # lower of LHS
    UdQ = zeros(cellxmax, cellymax, nval)               # upper of LHS

    RHS_temp = zeros(cellxmax, cellymax, nval)          # temporary RHS

    # define at cell boundaries
    lambda_facex = zeros(cellxmax+1, cellymax)          # lambda for computational time steps
    lambda_facey = zeros(cellxmax, cellymax+1)          # lambda for computational time steps

    # misc
    norm2 = zeros(nval)                                 # Residuals by norm-2
    I = zeros(nval, nval)                               # identity matrix
    for l in 1:nval
        I[l,l] = 1.0
    end

    return Qbasen, Qconn, Qconn_hat, Qbasem, dtau, lambda_facex, lambda_facey,
            A_adv_hat_m, A_adv_hat_p, B_adv_hat_m, B_adv_hat_p, A_beta_shig, B_beta_shig,
            jalphaP, jbetaP, delta_Q, delta_Q_temp, D, Lx, Ly, Ux, Uy, LdQ, UdQ, RHS_temp,
            norm2, I
end