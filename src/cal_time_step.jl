# ------------------------------------
# time integration with explicit methods
# ------------------------------------
function time_integration_explicit(dt, Qcon_hat, RHS, cellxmax, cellymax, nval)
    Threads.@threads for l in 1:nval
        for j in 2:cellymax-1
            for i in 2:cellxmax-1
                Qcon_hat[i,j,l] = Qcon_hat[i,j,l] + dt*RHS[i,j,l]
            end
        end
    end
    return Qcon_hat
end 

# ------------------------------------
# Calculation of the advection Jacobian matrix by one-wave approximation for lusgs
# ------------------------------------
function one_wave(A_adv_hat_m, A_adv_hat_p, B_adv_hat_m, B_adv_hat_p, A_beta_shig, B_beta_shig, I,
                    Qbase, Qcon, cellxmax, cellymax, vecAx, vecAy, specific_heat_ratio, volume, nval)
    beta = 1.1

    for j in 2:cellymax
        for i in 2:cellxmax
            for k in 1:2 #A,B            
                jacob_temp = zeros(nval,nval)
                if k == 1
                    kx_av = 0.5*(vecAx[i,j,1]+vecAx[i+1,j,1]) / volume[i,j]
                    ky_av = 0.5*(vecAx[i,j,2]+vecAx[i+1,j,2]) / volume[i,j]
                elseif k == 2
                    kx_av = 0.5*(vecAy[i,j,1]+vecAy[i,j+1,1]) / volume[i,j]
                    ky_av = 0.5*(vecAy[i,j,2]+vecAy[i,j+1,2]) / volume[i,j]
                end

                rho = Qbase[i,j,1]
                u = Qbase[i,j,2]
                v = Qbase[i,j,3]
                e = Qcon[i,j,4]
                p = Qbase[i,j,4]

                g = specific_heat_ratio
                
                Z = kx_av*u+ky_av*v
                q2   = 0.5*(u^2+v^2)
                b1c2 = 0.5*q2*(g-1)
                gebyrho = g*e/rho

                # advection Jacobian matrix in the general coordinate system
                jacob_temp[1,1] = 0.0
                jacob_temp[1,2] = kx_av
                jacob_temp[1,3] = ky_av
                jacob_temp[1,4] = 0.0
            
                jacob_temp[2,1] = -u*Z + kx_av*b1c2
                jacob_temp[2,2] = Z - (g-2)*kx_av*u
                jacob_temp[2,3] = ky_av*u - kx_av*(g-1)*v
                jacob_temp[2,4] = (g-1)*kx_av
            
                jacob_temp[3,1] = -v*Z + ky_av*b1c2
                jacob_temp[3,2] = kx_av*v - ky_av*(g-1)*u
                jacob_temp[3,3] = Z - ky_av*(g-2)*v
                jacob_temp[3,4] = (g-1)*ky_av
            
                jacob_temp[4,1] = Z*(-gebyrho + 2*b1c2)
                jacob_temp[4,2] = kx_av*(gebyrho-b1c2) - (g-1)*u*Z
                jacob_temp[4,3] = ky_av*(gebyrho-b1c2) - (g-1)*v*Z
                jacob_temp[4,4] = g*Z

                c = (g*rho/p)^0.5
                shigma = abs(Z) + c*(kx_av^2+ky_av^2)^0.5

                for l in 1:nval
                    I[l,l] = beta * shigma
                end

                if k == 1
                    for m in 1:nval
                        for l in 1:nval                        
                            A_adv_hat_p[i,j,l,m] = 0.5*(jacob_temp[l,m] + I[l,m])
                            A_adv_hat_m[i,j,l,m] = 0.5*(jacob_temp[l,m] - I[l,m])
                        end
                    end
                    A_beta_shig[i,j] = beta * shigma
                elseif k ==2                    
                    for m in 1:nval
                        for l in 1:nval
                            B_adv_hat_p[i,j,l,m] = 0.5*(jacob_temp[l,m] + I[l,m])
                            B_adv_hat_m[i,j,l,m] = 0.5*(jacob_temp[l,m] - I[l,m])
                        end
                    end
                    B_beta_shig[i,j] = beta * shigma
                end

            end
        end
    end
    
    # reset
    for l in 1:nval
        I[l,l] = 1.0
    end

    return A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, A_beta_shig, B_beta_shig
end

# ------------------------------------
# Calculation of the viscosity Jacobian matrix by center differential for lusgs
# ------------------------------------
function central_diff_jacobian(jalphaP, jbetaP, Qbase, Qcon, cellxmax, cellymax, mu, lambda,
                            vecAx, vecAy, specific_heat_ratio, volume, nval)
    
    # substitution to Corner cell
    Qbase[1,1,1] = 1.0
    Qbase[cellxmax,1,1] = 1.0
    Qbase[1,cellymax,1] = 1.0
    Qbase[cellxmax,cellymax,1] = 1.0
    volume[1,1,1] = 1.0
    volume[cellxmax,1,1] = 1.0
    volume[1,cellymax,1] = 1.0
    volume[cellxmax,cellymax,1] = 1.0
    
    for j in 1:cellymax
        for i in 1:cellxmax
            rho = Qbase[i,j,1]
            xi_x_av = 0.5*(vecAx[i,j,1] + vecAx[i+1,j,1]) / volume[i,j]
            xi_y_av = 0.5*(vecAx[i,j,2] + vecAx[i+1,j,2]) / volume[i,j]
            alpha = (xi_x_av^2+xi_y_av^2)^0.5
            jalphaP[i,j] = alpha / volume[i,j] * (2*mu[i,j]/rho)

            eta_x_av = 0.5*(vecAy[i,j,1] + vecAy[i,j+1,1]) / volume[i,j]
            eta_y_av = 0.5*(vecAy[i,j,2] + vecAy[i,j+1,2]) / volume[i,j]
            beta = (xi_x_av^2+xi_y_av^2)^0.5
            jbetaP[i,j] = beta / volume[i,j] * (2*mu[i,j]/rho)
        end
    end
    return jalphaP, jbetaP
end

# ------------------------------------
# lusgs
# ------------------------------------
function lusgs(D, Lx, Ly, Ux, Uy, LdQ, UdQ, RHS_temp, I, dt, dtau, Qcon_hat, Qconn_hat, delta_Q,
                A_adv_hat_p,  A_adv_hat_m,  B_adv_hat_p,  B_adv_hat_m,  A_beta_shig,  B_beta_shig,
                jalphaP,  jbetaP, RHS, cellxmax, cellymax, volume, nval)
       
    # calculate L and U
    for m in 1:nval
        for l in 1:nval
            for j in 1:cellymax
                for i in 1:cellxmax            
                    Lx[i,j,l,m] = dt*(A_adv_hat_p[i,j,l,m] + jalphaP[i,j]*I[l,m])
                    Ly[i,j,l,m] = dt*(B_adv_hat_p[i,j,l,m] + jbetaP[i,j]*I[l,m])
                    Ux[i,j,l,m] = dt*(A_adv_hat_m[i,j,l,m] - jalphaP[i,j]*I[l,m])
                    Uy[i,j,l,m] = dt*(B_adv_hat_m[i,j,l,m] - jbetaP[i,j]*I[l,m])
                end
            end
        end
    end

    for j in 2:cellymax-1
        for i in 2:cellxmax-1
            for l in 1:nval
                for m in 1:nval
                    LdQ[i,j,l] += Lx[i-1,j,l,m]*delta_Q[i-1,j,m] + Ly[i,j-1,l,m]*delta_Q[i,j-1,m]
                    UdQ[i,j,l] += Ux[i+1,j,l,m]*delta_Q[i+1,j,m] + Uy[i,j+1,l,m]*delta_Q[i,j+1,m]
                end
            end
        end
    end

    # diagonal
    for j in 2:cellymax-1
        for i in 2:cellxmax-1
            D[i,j] = dt/dtau[i,j] + 1.0 + dt*(A_beta_shig[i,j]+2*jalphaP[i,j] + B_beta_shig[i,j]+2*jbetaP[i,j])
        end
    end
    
    # RHS
    for l in 1:nval    
        for j in 2:cellymax-1
            for i in 2:cellxmax-1            
                RHS_temp[i,j,l] = - (Qcon_hat[i,j,l]-Qconn_hat[i,j,l])*1.0 + dt*RHS[i,j,l]
            end
        end
    end

    # lower sweep
    for l in 1:nval
        for j in 2:cellymax-1
            for i in 2:cellxmax-1
                delta_Q[i,j,l] = D[i,j]^(-1) * (RHS_temp[i,j,l]+LdQ[i,j,l])
            end
        end
    end             
    
    # upepr sweep
    for l in 1:nval
        for j in 2:cellymax-1
            for i in 2:cellxmax-1
                delta_Q[i,j,l] = delta_Q[i,j,l] - D[i,j]^(-1) * UdQ[i,j,l]
            end
        end
    end

    # reset
    for l in 1:nval
        for j in 2:cellymax-1
            for i in 2:cellxmax-1
                LdQ[i,j,l] = 0.0
                UdQ[i,j,l] = 0.0                
            end
        end
    end
    
    return delta_Q
end

