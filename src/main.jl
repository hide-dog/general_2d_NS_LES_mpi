using ProgressMeter
using Dates

function main()
    
    # Start of time measurement
    start_t = now()

    out_dir  = "result"         # output dir
    PARAMDAT = "PARAMDAT.json"  # directory
    fwrite   = "write"          # write file 
    
    nval = 4         # number of conserved variables
    Rd   = 287.0     # gas constant of air, J/(kg K)
    R    = 8.314     # gas constant, J/(K mol)
    
    # read grids and parameter
    xmax, ymax, nodes, vecAx_all, vecAy_all = read_allgrid()
    out_file_front, out_ext, restartnum, restart_file, init_small, norm_ok,
    time_integ, nt, dt, every_outnum, in_nt, dtau, cfl, ad_scheme,
    init_rho, init_u, init_v, init_p, init_T, specific_heat_ratio, Rd, bdcon = input_para(PARAMDAT)
    
    # number of all cells
    cellxmax = xmax - 1
    cellymax = ymax - 1
    
    # allocation for all cells
    Qbase_all, volume_all, cellcenter_all, wally_all, dx_all, dy_all = common_mpi_allocation(cellxmax, cellymax, nval)
    
    # set initial condition
    Qbase_all, restartnum = set_initQbase(Qbase_all, cellxmax, cellymax, restart_file, init_rho, init_u, init_v, init_p, init_T,
                                      specific_heat_ratio, out_file_front, out_ext, out_dir, restartnum, Rd, nval)
    
    # set volume, dx and dy
    volume_all     = set_volume(nodes, cellxmax, cellymax, volume_all)
    cellcenter_all = set_cellcenter(cellcenter_all, nodes, cellxmax, cellymax)
    dx_all, dy_all = set_dx_lts(dx_all, dy_all, nodes, cellxmax, cellymax)
    reset_write(fwrite)

    # wally
    wally_all, swith_wall = set_wally(nodes, bdcon, wally_all, cellcenter_all, cellxmax, cellymax)

    # AUSM+up
    Minf = 0.0
    if ad_scheme == 2
        Minf = set_Minf(bdcon, specific_heat_ratio, Rd, nval)
    end

    # mpi set up
    MPI.Init()
    cellNx, ista, iend, cellNy, jsta, jend, num_div, rank, size, comm = set_mpi(cellxmax, cellymax)
    
    # number of mpi cells + virtual cell
    cellx_mpi = cellNx + 4
    celly_mpi = cellNy + 4
    
    # mpi_allocation
    Qbase, volume, cellcenter, wally, yplus, dx, dy, Qcon, Qcon_hat, mu, lambda, 
    vecAx, vecAy, E_adv_hat, F_adv_hat, E_vis_hat, F_vis_hat, RHS    = common_mpi_allocation(cellx_mpi, celly_mpi, nval)

    Qbase, volume, cellcenter, wally, dx, dy = set_value_for_mpi(Qbase, cellx_mpi, celly_mpi, volume, cellcenter, wally, dx, dy,
                                                    Qbase_all, volume_all, cellcenter_all, wally_all, dx_all, dy_all, nval, num_div, rank)
    
    vecAx, vecAy = set_vecA_mpi(vecAx_all, vecAy_all, vecAx, vecAy, cellx_mpi, celly_mpi, num_div, rank)

    boundary_number = set_boundary_number(bdcon, rank, size)
    
    # check boundary condition
    check_bd(bdcon)

    # set initial condition for imaginary cell
    Qbase    = set_boundary(Qbase, cellx_mpi, celly_mpi, vecAx, vecAy, bdcon, Rd, specific_heat_ratio, 
                            nval, boundary_number, rank, size, comm)

    # write number of threads
    print("threads num : ")
    println(Threads.nthreads())

    # main loop
    loop_ite = 0
    if time_integ == "1"
        # exlicit scheme
        prog = Progress(nt,1)
        @time for t in 1:nt
            next!(prog)

            loop_ite += 1
            
            # step number
            evalnum = t + restartnum
            
            # set conserved variables in the general coordinate system
            Qbase    = set_boundary(Qbase, cellx_mpi, celly_mpi, vecAx, vecAy, bdcon, Rd, specific_heat_ratio, 
                                    nval, boundary_number, rank, size, comm)
            Qcon     = base_to_conservative(Qbase, Qcon, cellxmax, cellymax, specific_heat_ratio)
            Qcon_hat = setup_Qcon_hat(Qcon, Qcon_hat, cellxmax, cellymax, volume, nval)
            
            # set viscosity and thermal Conductivity
            mu     = set_mu(mu, Qbase, cellxmax, cellymax, specific_heat_ratio, Rd)
            lambda = set_lambda(lambda, Qbase, cellxmax, cellymax, mu, specific_heat_ratio, Rd)
            
            # yplus
            yplus = cal_yplus(yplus, Qbase, wally, swith_wall, mu, cellxmax, cellymax, vecAx, vecAy, volume)
            #yplus = ones(cellxmax, cellymax)*100                  # yplus
            
            # advection_term
            E_adv_hat, F_adv_hat = AUSM(E_adv_hat, F_adv_hat, Qbase, Qcon, cellxmax, cellymax, vecAx, vecAy, 
                                        specific_heat_ratio, volume, nval, Minf, ad_scheme)
                        
            # viscos_term
            E_vis_hat, F_vis_hat = central_diff(E_vis_hat, F_vis_hat, Qbase, Qcon, cellxmax, cellymax, mu, lambda,
                                                vecAx, vecAy, specific_heat_ratio, volume, Rd, nval, yplus, swith_wall)
            
            #println(" fff ")
            #println(E_adv_hat[10,10,:])
            #println(F_adv_hat[10,10,:])
            #println(E_vis_hat[10,10,:])
            #println(F_vis_hat[10,10,:])
            

            #throw(UndefVarError(:x))
            
            #println(yplus[:,1])
            #println(" yplus ")
            #println(yplus[1000,2])
            #println(yplus[1000,3])
            #println(mu[1000,2] / Qbase[1000,2,1])
            #println(yplus[:,4])
            #println(yplus[:,5])
            #println(wally[1000,2])
            #println(wally[1000,3])
            #throw(UndefVarError(:x))
            

            # RHS
            RHS = setup_RHS(RHS, cellxmax, cellymax, E_adv_hat, F_adv_hat, E_vis_hat, F_vis_hat, nval, volume)
            
            # time integral
            Qcon_hat = time_integration_explicit(dt, Qcon_hat, RHS, cellxmax, cellymax, nval)
            
            # calculate primitive variables
            Qcon  = Qhat_to_Q(Qcon, Qcon_hat, cellxmax, cellymax, volume, nval)
            Qbase = conservative_to_base(Qbase, Qcon, cellxmax, cellymax, specific_heat_ratio)
            #Qbase_ave = cal_Qave(Qbase, Qbase_ave, cellxmax, cellymax, nval)
            
            # output
            if round(evalnum) % every_outnum == 0
                println("\n")
                println("nt_______________________________"*string(round(evalnum)))
                output_result(evalnum, Qbase, cellxmax, cellymax, specific_heat_ratio, out_file_front, out_ext, out_dir, Rd, nval)
                #output_ave(Qbase_ave, cellxmax, cellymax, out_file_front, out_ext, out_dir, Rd, nval, loop_ite)
            end
            
            # Find out if the results were divergent
            check_divrege(Qbase, cellxmax, cellymax, Rd, fwrite)
        end
    elseif time_integ == "2"
        # implicit

        # allocation for implicit
        Qbasen, Qconn, Qconn_hat, Qbasem, dtau, lambda_facex, lambda_facey,
        A_adv_hat_m, A_adv_hat_p, B_adv_hat_m, B_adv_hat_p, A_beta_shig, B_beta_shig,
        jalphaP, jbetaP, delta_Q, delta_Q_temp, D, Lx, Ly, Ux, Uy, LdQ, UdQ, RHS_temp,
        norm2, I = allocation_implicit(cellxmax, cellymax, nval)

        prog = Progress(nt,1)
        @time for t in 1:nt
            next!(prog)
            
            # step number
            evalnum = t + restartnum
            
            # write physical time 
            output_physicaltime(fwrite, t, dt)
            
            # copy
            for l in 1:nval
                for j in 1:cellymax
                    for i in 1:cellxmax
                        Qbasen[i,j,l] = Qbase[i,j,l]
                        Qbasem[i,j,l] = Qbase[i,j,l]
                    end
                end
            end

            # set conserved variables in the general coordinate system for inner iteration
            Qconn     = base_to_conservative(Qbasen, Qconn, cellxmax, cellymax, specific_heat_ratio)
            Qconn_hat = setup_Qcon_hat(Qconn, Qconn_hat, cellxmax, cellymax, volume, nval)     

            # start inner iteration
            for tau in 1:in_nt
                # set conserved variables in the general coordinate system
                Qbasem = set_boundary(Qbasem, cellxmax, cellymax, vecAx, vecAy, bdcon, Rd, specific_heat_ratio, nval)
                Qcon     = base_to_conservative(Qbasem, Qcon, cellxmax, cellymax, specific_heat_ratio)
                Qcon_hat = setup_Qcon_hat(Qcon, Qcon_hat, cellxmax, cellymax, volume, nval)

                # set viscosity and thermal Conductivity
                mu     = set_mu(mu, Qbasem, cellxmax, cellymax, specific_heat_ratio, Rd)
                lambda = set_lambda(lambda, Qbasem, cellxmax, cellymax, mu, specific_heat_ratio, Rd)
                
                yplus = cal_yplus(yplus, Qbasem, wally, swith_wall, mu, cellxmax, cellymax, vecAx, vecAy, volume)
                #println(yplus[:,2])
            
                #throw(UndefVarError(:x))

                # set inner time step by local time stepping
                dtau   = set_lts(dtau, lambda_facex, lambda_facey, Qbasem, cellxmax, cellymax, mu, dx, dy,
                                vecAx, vecAy, volume, specific_heat_ratio, cfl)
                                
                #advection_term
                E_adv_hat, F_adv_hat = AUSM(E_adv_hat, F_adv_hat, Qbasem, Qcon, cellxmax, cellymax, vecAx, vecAy, 
                                            specific_heat_ratio, volume, nval, Minf, ad_scheme)

                # viscos_term
                E_vis_hat, F_vis_hat = central_diff(E_vis_hat, F_vis_hat, Qbasem, Qcon, cellxmax, cellymax, mu, lambda,
                                                vecAx, vecAy, specific_heat_ratio, volume, Rd, nval, yplus, swith_wall)
                
                # RHS
                RHS = setup_RHS(RHS, cellxmax, cellymax, E_adv_hat, F_adv_hat, E_vis_hat, F_vis_hat, nval, volume)
            
                # lusgs_advection_term
                A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, 
                A_beta_shig, B_beta_shig = one_wave(A_adv_hat_m, A_adv_hat_p, B_adv_hat_m, B_adv_hat_p, A_beta_shig, B_beta_shig, I,
                                                    Qbase, Qcon, cellxmax, cellymax, vecAx, vecAy, specific_heat_ratio, volume, nval)
                # lusgs_viscos_term
                jalphaP, jbetaP = central_diff_jacobian(jalphaP, jbetaP, Qbase, Qcon, cellxmax, cellymax, mu, lambda,
                                                        vecAx, vecAy, specific_heat_ratio, volume, nval)

                #println(RHS[2,:,1])
                #println(A_adv_hat_m[2,:,1])
                #println(jalphaP[2,:])
                
                # LUSGS
                ite = 0
                while true
                    # copy delta_Q
                    for l in 1:nval
                        for j in 1:cellymax
                            for i in 1:cellxmax
                                delta_Q_temp[i,j,l] = delta_Q[i,j,l]
                            end
                        end
                    end
                                        
                    # Reversing the left-hand side by lusgs
                    delta_Q = lusgs(D, Lx, Ly, Ux, Uy, LdQ, UdQ, RHS_temp, I, dt, dtau, Qcon_hat, Qconn_hat, delta_Q,
                                    A_adv_hat_p,  A_adv_hat_m,  B_adv_hat_p,  B_adv_hat_m,  A_beta_shig,  B_beta_shig,
                                     jalphaP,  jbetaP, RHS, cellxmax, cellymax, volume, nval)
                    
                    # cal Residuals by norm-2
                    res   = set_res(delta_Q, delta_Q_temp, cellxmax, cellymax, nval)
                    norm2 = check_converge(res, RHS, cellxmax, cellymax, init_small, nval)
                    
                    # Find out if the results converged
                    if norm2[1] < norm_ok && norm2[2] < norm_ok && norm2[3] < norm_ok && norm2[4] < norm_ok
                        break
                    end
                    if ite % 100 ==0
                        println(" now cal norm2 ")
                        println(norm2)
                        if isequal(norm2[1], NaN) == true
                            println("  ")
                            println(" norm2 = NaN ")
                            println(" stop cal ")
                            throw(UndefVarError(:x))
                        end
                    end

                    ite += 1
                end

                # output inner time
                output_innertime(fwrite, tau, norm2, nval)
                
                # Updating
                for l in 1:nval
                    for j in 2:cellymax-1
                        for i in 2:cellxmax-1
                            Qcon_hat[i,j,l] = Qcon_hat[i,j,l] + delta_Q[i,j,l]
                        end
                    end
                end
                
                # calculate primitive variables
                Qcon = Qhat_to_Q(Qcon, Qcon_hat, cellxmax, cellymax, volume, nval)
                Qbasem = conservative_to_base(Qbasem, Qcon, cellxmax, cellymax, specific_heat_ratio)

                # Find out if the results were divergent
                check_divrege(Qbasem, cellxmax, cellymax, Rd, fwrite)
            end
            # End of the inner iteration

            # Updating in physical time
            for l in 1:nval
                for j in 1:cellymax
                    for i in 1:cellxmax
                        Qbase[i,j,l] = Qbasem[i,j,l]
                    end
                end
            end
            
            # output
            if round(evalnum) % every_outnum == 0
                println("\n")
                println("nt_______________________________"*string(round(evalnum)))
                output_result(evalnum, Qbase, cellxmax, cellymax, specific_heat_ratio, out_file_front, out_ext, out_dir, Rd, nval)
            end

            # Find out if the results were divergent
            check_divrege(Qbase, cellxmax, cellymax, Rd, fwrite)
        end
    end
    
    # end of mpi
    MPI.Barrier(comm)
    MPI.Finalize() 

    # end of time measurement
    end_t = now()

    # output of calculation time
    output_fin(fwrite, start_t, end_t, nt, dt, in_nt, cellxmax, cellymax)
end


# -- main --
main()
#throw(UndefVarError(:x))