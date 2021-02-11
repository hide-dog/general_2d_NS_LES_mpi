# ------------------------------------
# calcucation of the viscosity term by central difference
# ------------------------------------
function central_diff(E_vis_hat, F_vis_hat, Qbase, Qcon, cellxmax, cellymax, mu, lambda,
                        vecAx, vecAy, specific_heat_ratio, volume, Rd, nval, yplus, swith_wall)
    
    for j in 2:cellymax -1
        for i in 2:cellxmax+1 -1
            rho_av = 0.5*(Qbase[i,j,1] + Qbase[i-1,j,1]) 
            u_av   = 0.5*(Qbase[i,j,2] + Qbase[i-1,j,2])
            v_av   = 0.5*(Qbase[i,j,3] + Qbase[i-1,j,3])
            
            mu_av     = 0.5*(mu[i-1,j] + mu[i,j]) 
            lambda_av = (0.5*(1/lambda[i-1,j] + 1/lambda[i,j]))^(-1)

            volume_av = 0.5*(volume[i,j] + volume[i-1,j])
    
            dudxi = Qbase[i,j,2] - Qbase[i-1,j,2]
            dvdxi = Qbase[i,j,3] - Qbase[i-1,j,3]
            dTdxi = Qbase[i,j,4]/(Qbase[i,j,1]*Rd) - Qbase[i-1,j,4]/(Qbase[i-1,j,1]*Rd)
    
            dudeta = 0.25 * (Qbase[i,j+1,2] - Qbase[i,j-1,2] + Qbase[i-1,j+1,2] - Qbase[i-1,j-1,2])
            dvdeta = 0.25 * (Qbase[i,j+1,3] - Qbase[i,j-1,3] + Qbase[i-1,j+1,3] - Qbase[i-1,j-1,3])
            dTdeta = 0.25 * (Qbase[i,j+1,4]/(Qbase[i,j+1,1]*Rd) - Qbase[i,j-1,4]/(Qbase[i,j-1,1]*Rd) +
                            Qbase[i-1,j+1,4]/(Qbase[i-1,j+1,1]*Rd) -  Qbase[i-1,j-1,4]/(Qbase[i-1,j-1,1]*Rd))           

            vecAy_xav    = 0.25*( vecAy[i,j,1] + vecAy[i,j+1,1] + vecAy[i-1,j,1] + vecAy[i-1,j+1,1] )
            vecAy_yav    = 0.25*( vecAy[i,j,2] + vecAy[i,j+1,2] + vecAy[i-1,j,2] + vecAy[i-1,j+1,2] )

            dudx = vecAx[i,j,1] * dudxi + vecAy_xav * dudeta
            dvdy = vecAx[i,j,2] * dvdxi + vecAy_yav * dvdeta
            dudy = vecAx[i,j,2] * dudxi + vecAy_yav * dudeta
            dvdx = vecAx[i,j,1] * dvdxi + vecAy_xav * dvdeta
            
            sigma_xx = 2/3*mu_av*( 2*dudx - dvdy )          
                    
            sigma_yy = 2/3*mu_av*( -dudx  + 2*dvdy )          
                    
            sigma_xy = mu_av*( dudy + dvdx )
            
            dTdx = vecAx[i,j,1]*dTdxi + vecAy_xav*dTdeta
            dTdy = vecAx[i,j,2]*dTdxi + vecAy_yav*dTdeta
       
            betax = u_av*sigma_xx + v_av*sigma_xy + lambda_av * dTdx
            betay = u_av*sigma_xy + v_av*sigma_yy + lambda_av * dTdy
            
            yplus_av = 0.5 *(yplus[i-1,j] + yplus[i,j])

            
            if swith_wall[1] == 1 && i == 2
                tau_xx, tau_xy, tau_yy, e_sgs_x, e_sgs_y = 0, 0, 0, 0, 0
            elseif swith_wall[2] == 1 && i == cellxmax
                tau_xx, tau_xy, tau_yy, e_sgs_x, e_sgs_y = 0, 0, 0, 0, 0
            else
                tau_xx, tau_xy, tau_yy, e_sgs_x, e_sgs_y = Smagorinsky_model(dudx, dvdy, dudy, dvdx, dTdx, dTdy, rho_av, u_av, v_av, volume_av, yplus_av,i,j)
                #tau_xx, tau_xy, tau_yy, e_sgs_x, e_sgs_y = 0, 0, 0, 0, 0

            #=    
                if i == 51 && j ==3
                    println(" tur tau ")
                    println(tau_xx)
                    println(tau_xy)
                    println(tau_yy)
                    println(e_sgs_x)
                    println(e_sgs_y)
                    println(" kiyo ")
                    println(vecAx[i,j,1]*sigma_xx + vecAx[i,j,2]*sigma_xy)
                    println(vecAx[i,j,1]*tau_xx + vecAx[i,j,2]*tau_xy)
                    println(vecAx[i,j,1]*tau_xx)
                    println(vecAx[i,j,2]*tau_xy)
                end
                =#

            end
            
            E_vis_hat[i,j,1] = 0.0
            E_vis_hat[i,j,2] = ((vecAx[i,j,1]*sigma_xx + vecAx[i,j,2]*sigma_xy) - (vecAx[i,j,1]*tau_xx + vecAx[i,j,2]*tau_xy)) / volume_av
            E_vis_hat[i,j,3] = ((vecAx[i,j,1]*sigma_xy + vecAx[i,j,2]*sigma_yy) - (vecAx[i,j,1]*tau_xy + vecAx[i,j,2]*tau_yy)) / volume_av
            E_vis_hat[i,j,4] = ((vecAx[i,j,1]*betax + vecAx[i,j,2]*betay) - (vecAx[i,j,1]*e_sgs_x + vecAx[i,j,2]*e_sgs_y))   / volume_av
        end
    end
    
    for j in 2:cellymax+1 -1
        for i in 2:cellxmax -1
            rho_av = 0.5*(Qbase[i,j,1] + Qbase[i,j-1,1]) 
            u_av = 0.5*(Qbase[i,j,2] + Qbase[i,j-1,2])
            v_av = 0.5*(Qbase[i,j,3] + Qbase[i,j-1,3])

            mu_av     = 0.5*(mu[i,j-1] + mu[i,j])
            lambda_av = (0.5*(1/lambda[i,j-1] + 1/lambda[i,j]))^(-1)

            volume_av = 0.5*(volume[i,j] + volume[i,j-1])
    
            dudxi = 0.25 * (Qbase[i+1,j,2] - Qbase[i-1,j,2] + Qbase[i+1,j-1,2] - Qbase[i-1,j-1,2])
            dvdxi = 0.25 * (Qbase[i+1,j,3] - Qbase[i-1,j,3] + Qbase[i+1,j-1,3] - Qbase[i-1,j-1,3])
            dTdxi = 0.25 * (Qbase[i+1,j,4]/(Qbase[i+1,j,1]*Rd) - Qbase[i-1,j,4]/(Qbase[i-1,j,1]*Rd) +
                            Qbase[i+1,j-1,4]/(Qbase[i+1,j-1,1]*Rd) - Qbase[i-1,j-1,4]/(Qbase[i-1,j-1,1]*Rd))

            dudeta = Qbase[i,j,2] - Qbase[i,j-1,2]
            dvdeta = Qbase[i,j,3] - Qbase[i,j-1,3]
            dTdeta = (Qbase[i,j,4]/(Qbase[i,j,1]*Rd) - Qbase[i,j-1,4]/(Qbase[i,j-1,1]*Rd))
        
            vecAx_xav    = 0.25*( vecAx[i,j,1] + vecAx[i+1,j,1] + vecAx[i,j-1,1] + vecAx[i+1,j-1,1] )
            vecAx_yav    = 0.25*( vecAx[i,j,2] + vecAx[i+1,j,2] + vecAx[i,j-1,2] + vecAx[i+1,j-1,2] )

            dudx = vecAx_xav * dudxi + vecAy[i,j,1] * dudeta
            dvdy = vecAx_yav * dvdxi + vecAy[i,j,2] * dvdeta
            dudy = vecAx_yav * dudxi + vecAy[i,j,2] * dudeta
            dvdx = vecAx_xav * dvdxi + vecAy[i,j,1] * dvdeta
        
            sigma_xx = 2/3*mu_av*(2*dudx - dvdy)
                    
            sigma_yy = 2/3*mu_av*(-dudx + 2*dvdy)
                    
            sigma_xy = mu_av*(dudy + dvdx)
            
            dTdx = vecAx_xav*dTdxi + vecAy[i,j,1]*dTdeta
            dTdy = vecAx_yav*dTdxi + vecAy[i,j,2]*dTdeta
                   
            betax = u_av*sigma_xx + v_av*sigma_xy + lambda_av * dTdx
            betay = u_av*sigma_xy + v_av*sigma_yy + lambda_av * dTdy
            
            yplus_av = 0.5 *(yplus[i,j-1] + yplus[i,j])

            if swith_wall[3] == 1 && j == 2
                tau_xx, tau_xy, tau_yy, e_sgs_x, e_sgs_y = 0, 0, 0, 0, 0
            elseif swith_wall[4] == 1 && j == cellymax
                tau_xx, tau_xy, tau_yy, e_sgs_x, e_sgs_y = 0, 0, 0, 0, 0
            else
                #tau_xx, tau_xy, tau_yy, e_sgs_x, e_sgs_y = Smagorinsky_model(dudx, dvdy, dudy, dvdx, dTdx, dTdy, rho_av, u_av, v_av, volume_av, yplus_av,i,j)
                tau_xx, tau_xy, tau_yy, e_sgs_x, e_sgs_y = 0, 0, 0, 0, 0
            end
            
            F_vis_hat[i,j,1] = 0.0
            F_vis_hat[i,j,2] = ((vecAy[i,j,1]*sigma_xx + vecAy[i,j,2]*sigma_xy) - (vecAy[i,j,1]*tau_xx + vecAy[i,j,2]*tau_xy)) / volume_av
            F_vis_hat[i,j,3] = ((vecAy[i,j,1]*sigma_xy + vecAy[i,j,2]*sigma_yy) - (vecAy[i,j,1]*tau_xy + vecAy[i,j,2]*tau_yy)) / volume_av
            F_vis_hat[i,j,4] = ((vecAy[i,j,1]*betax + vecAy[i,j,2]*betay) - (vecAy[i,j,1]*e_sgs_x + vecAy[i,j,2]*e_sgs_y))  / volume_av
        end
    end
    return E_vis_hat, F_vis_hat
end

function Smagorinsky_model(dudx, dvdy, dudy, dvdx, dTdx, dTdy, rho_av, u_av, v_av, volume_av, yplus_av,i,j)
    Cs = 0.18 
    Cl = 0.09
    Pr_sgs = 0.6
    k = 0.41

    S11 = dudx
    S12 = 0.5*(dudy + dvdx)
    S21 = 0.5*(dudy + dvdx)
    S22 = dvdy

    absS = ( 2 * (S11^2 + S12^2 + S21^2 + S22^2) )^0.5
    Delta = (volume_av) ^ (1/3) * wallf_Van_Driest(yplus_av)

    #lsgs = 0.5*((Cs*Delta - k*y) - abs(Cs*Delta - k*y)) + k*y # min([ky, Cs*Delta])
    #nu_sgs = (lsgs)^2 * absS
    nu_sgs = (Cs*Delta)^2 * absS

    # SGS stress
    #tau_xx = -2 * rho_av * nu_sgs * 2/3 * dudx +  1/3 * Cl * 2 * rho_av * Delta^2 * absS^2
    #tau_xy = - rho_av * nu_sgs * (dvdx + dudy)
    #tau_yy = -2 * rho_av * nu_sgs * 2/3 * dvdy +  1/3 * Cl * 2 * rho_av * Delta^2 * absS^2
    tau_xx = -2 * rho_av * nu_sgs * S11
    tau_xy = -2 * rho_av * nu_sgs * S12
    tau_yy = -2 * rho_av * nu_sgs * S22

    # SGS energy
    e_sgs_x = -rho_av * nu_sgs / Pr_sgs * dTdx + tau_xx*u_av + tau_xy*v_av
    e_sgs_y = -rho_av * nu_sgs / Pr_sgs * dTdy + tau_xy*u_av + tau_yy*v_av

    # 
    tau_xx = tau_xx / volume_av   *1.0
    tau_xy = tau_xy / volume_av   *1.0
    tau_yy = tau_yy / volume_av   *1.0
    e_sgs_x = e_sgs_x / volume_av *1.0
    e_sgs_y = e_sgs_y / volume_av *1.0
    
    #=
    if i == 51 && j ==3
        println(" iroiro ")
        println(absS)
        println(nu_sgs)
        println(-2 * rho_av * nu_sgs * 2/3 * dudx)
        println(1/3 * Cl * 2 * rho_av * Delta^2 * absS^2)
    end
    =#

    return tau_xx, tau_xy, tau_yy, e_sgs_x, e_sgs_y
end