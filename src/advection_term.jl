# ------------------------------------
# advection term by AUSM
# ------------------------------------
function AUSM(E_adv_hat, F_adv_hat, Qbase, Qcon, cellxmax, cellymax, vecAx, vecAy, specific_heat_ratio, volume, nval, Minf, ad_scheme)
    g         = specific_heat_ratio
    temp_vecX = zeros(nval)   # vecAx used for the pressure term
    temp_vecY = zeros(nval)   # vecAx used for the pressure term
    Lpsi = zeros(nval)        
    Rpsi = zeros(nval)

    mdot = 0.0
    ph = 0.0
    
    tvAxx = 0.0
    tvAxy = 0.0
    vAxx = 0.0
    vAxy = 0.0
    tvAyx = 0.0
    tvAyy = 0.0
    vAyx = 0.0
    vAyy = 0.0
    for j in 2:cellymax -1
        for i in 2:cellxmax+1 -1
            volume_av = 0.5*(volume[i-1,j] + volume[i,j])

            # i-1 cell
            # 規格化
            tvAxx = 0.5*(vecAx[i-1,j,1]+vecAx[i,j,1])
            tvAxy = 0.5*(vecAx[i-1,j,2]+vecAx[i,j,2])
            vAxx = tvAxx / (tvAxx^2+tvAxy^2)^0.5
            vAxy = tvAxy / (tvAxx^2+tvAxy^2)^0.5
            
            rhoL = Qbase[i-1,j,1]
            UL = Qbase[i-1,j,2]*vAxx + Qbase[i-1,j,3]*vAxy
            pL = Qbase[i-1,j,4]
            
            # i cell
            # 規格化
            tvAxx = 0.5*(vecAx[i,j,1]+vecAx[i+1,j,1])
            tvAxy = 0.5*(vecAx[i,j,2]+vecAx[i+1,j,2])
            vAxx = tvAxx / (tvAxx^2+tvAxy^2)^0.5
            vAxy = tvAxy / (tvAxx^2+tvAxy^2)^0.5
            
            rhoR = Qbase[i,j,1]
            UR = Qbase[i,j,2]*vAxx + Qbase[i,j,3]*vAxy
            pR = Qbase[i,j,4]
            
            #=
            if i==100 && j== 200
                mdot, ph = AUSM_plus_half(rhoL, rhoR, UL, UR, pL, pR, g,i,j)
                println("AUSM+")
                println(mdot)
                println(ph)
                
                mdot, ph = AUSM_plusup_half(rhoL, rhoR, UL, UR, pL, pR, Minf, g,i,j)
                println("AUSM+up")
                println(mdot)
                println(ph)
                
                velocity = 0.5 * ((Qbase[i-1,j,2]^2 + Qbase[i-1,j,3]^2)^0.5 + (Qbase[i,j,2]^2 + Qbase[i,j,3]^2)^0.5)
                mdot, ph = SLAU_half(rhoL, rhoR, UL, UR, pL, pR, velocity, specific_heat_ratio)
                println("SLAU")
                println(mdot)
                println(ph)
            end
            =#
        
            
            
            # scheme
            if  ad_scheme == 1
                mdot, ph = AUSM_plus_half(rhoL, rhoR, UL, UR, pL, pR, g,i,j)
            elseif ad_scheme == 2
                mdot, ph = AUSM_plusup_half(rhoL, rhoR, UL, UR, pL, pR, Minf, g, i,j)
            elseif ad_scheme == 4
                velocity = (0.5 * (Qbase[i-1,j,2]^2 + Qbase[i-1,j,3]^2 + Qbase[i,j,2]^2 + Qbase[i,j,3]^2))^0.5
                mdot, ph = SLAU_half(rhoL, rhoR, UL, UR, pL, pR, velocity, specific_heat_ratio)
            end
            
            # flux half
            sqAx = (vecAx[i,j,1]^2 + vecAx[i,j,2]^2)^0.5
            temp_vecX[2] = vecAx[i,j,1] / sqAx
            temp_vecX[3] = vecAx[i,j,2] / sqAx

            for l in 1:nval
                Lpsi[l] = Qcon[i-1,j,l] / Qcon[i-1,j,1]
                Rpsi[l] = Qcon[i,j,l]   / Qcon[i,j,1]
            end
            Lpsi[4] = (Qcon[i-1,j,4] + Qbase[i-1,j,4]) / Qcon[i-1,j,1]
            Rpsi[4] = (Qcon[i,j,4]   + Qbase[i,j,4])   / Qcon[i,j,1]

            if mdot > 0
                for l in 1:nval
                    E_adv_hat[i,j,l] = (mdot * Lpsi[l] + ph * temp_vecX[l]) * sqAx
                end
            else
                for l in 1:nval
                    E_adv_hat[i,j,l] = (mdot * Rpsi[l] + ph * temp_vecX[l]) * sqAx
                end
            end            
        end
    end

    for j in 2:cellymax+1 -1
        for i in 2:cellxmax -1
            volume_av = 0.5*(volume[i,j-1] + volume[i,j])

            # j-1 cell
            # 規格化
            tvAyx = 0.5*(vecAy[i,j-1,1]+vecAy[i,j,1])
            tvAyy = 0.5*(vecAy[i,j-1,2]+vecAy[i,j,2])
            vAyx = tvAyx / (tvAyx^2+tvAyy^2)^0.5
            vAyy = tvAyy / (tvAyx^2+tvAyy^2)^0.5
            
            rhoL = Qbase[i,j-1,1]
            VL = Qbase[i,j-1,2]*vAyx + Qbase[i,j-1,3]*vAyy
            pL = Qbase[i,j-1,4]
            
            # j cell
            tvAyx = 0.5*(vecAy[i,j,1]+vecAy[i,j+1,1])
            tvAyy = 0.5*(vecAy[i,j,2]+vecAy[i,j+1,2])
            vAyx = tvAyx / (tvAyx^2+tvAyy^2)^0.5
            vAyy = tvAyy / (tvAyx^2+tvAyy^2)^0.5
            
            rhoR = Qbase[i,j,1]
            VR = Qbase[i,j,2]*vAyx + Qbase[i,j,3]*vAyy
            pR = Qbase[i,j,4]

            # scheme
            if  ad_scheme == 1
                mdot, ph = AUSM_plus_half(rhoL, rhoR, VL, VR, pL, pR, g,i,j)
            elseif ad_scheme == 2
                mdot, ph = AUSM_plusup_half(rhoL, rhoR, VL, VR, pL, pR, Minf, g, i,j)
            elseif ad_scheme == 4
                velocity = (0.5 * (Qbase[i,j-1,2]^2 + Qbase[i,j-1,3]^2 + Qbase[i,j,2]^2 + Qbase[i,j,3]^2))^0.5
                mdot, ph = SLAU_half(rhoL, rhoR, VL, VR, pL, pR, velocity, specific_heat_ratio)
            end
            #=
            if i==170 && j== 2
                mdot, ph = AUSM_plus_half(rhoL, rhoR, VL,VR, pL, pR, g)
                println("AUSM+")
                println(mdot)
                println(ph)
                
                mdot, ph = AUSM_plusup_half(rhoL, rhoR, VL,VR, pL, pR, Minf, g)
                println("AUSM+up")
                println(mdot)
                println(ph)
                
                velocity = (0.5 * (Qbase[i,j-1,2]^2 + Qbase[i,j-1,3]^2 + Qbase[i,j,2]^2 + Qbase[i,j,3]^2))^0.5
                mdot, ph = SLAU_half(rhoL, rhoR, VL, VR, pL, pR, velocity, specific_heat_ratio)
                println("SLAU")
                println(mdot)
                println(ph)
            end
            =#
            
            # flux half
            sqAy = (vecAy[i,j,1]^2 + vecAy[i,j,2]^2)^0.5
            temp_vecY[2] = vecAy[i,j,1] / sqAy
            temp_vecY[3] = vecAy[i,j,2] / sqAy

            for l in 1:nval
                Lpsi[l] = Qcon[i,j-1,l] / Qcon[i,j-1,1]
                Rpsi[l] = Qcon[i,j,l]   / Qcon[i,j,1]
            end
            Lpsi[4] = (Qcon[i,j-1,4] + Qbase[i,j-1,4]) / Qcon[i,j-1,1]
            Rpsi[4] = (Qcon[i,j,4]   + Qbase[i,j,4])   / Qcon[i,j,1]
            
            if mdot > 0
                for l in 1:nval
                    F_adv_hat[i,j,l] = (mdot * Lpsi[l] + ph * temp_vecY[l]) * sqAy
                end
            else
                for l in 1:nval
                    F_adv_hat[i,j,l] = (mdot * Rpsi[l] + ph * temp_vecY[l]) * sqAy
                end
            end
        end
    end
    
    return E_adv_hat, F_adv_hat
end

function AUSM_plus_half(rhoL, rhoR, UL, UR, pL, pR, g,i,j)
    # param
    beta  = 1/8
    alpha = 3/16
    
    # L, R
    aL = (g * pL / rhoL)^0.5
    aR = (g * pR / rhoR)^0.5

    # half
    ah   = 0.5*( aL + aR )
    rhoh = 0.5*( rhoL + rhoR )

    ML = UL/ah
    MR = UR/ah

    M_p4 = 0
    M_m4 = 0
    p_p5 = 0
    p_m5 = 0
    if abs(ML) >= 1
        M_p4 = 0.5*(ML + abs(ML))
        p_p5 = 0.5*(ML + abs(ML)) / ML
    else
        M_p4 = 0.25*(ML+1)^2 + beta*(ML^2-1)^2
        p_p5 = 0.25*(ML+1)^2 * (2-ML) + alpha*ML*(ML^2-1)^2
    end
    
    if abs(MR) >= 1
        M_m4 = 0.5*(MR - abs(MR))
        p_m5 = 0.5*(MR - abs(MR)) / MR
    else
        M_m4 = -0.25*(MR-1)^2 - beta*(MR^2-1)^2
        p_m5 = 0.25*(MR-1)^2 * (2+MR) - alpha*MR*(MR^2-1)^2
    end
    
    # Mh = M_p4 + M_m4 + Mp/fa
    Mh = M_p4 + M_m4
#=
    if i==100 && j==200
        println(" MMM_AUSM+ ")
        println(M_p4)
        println(M_m4)
    end
  =#  
    # mdot half
    mdot = ah * Mh
    if Mh > 0
        mdot = mdot * rhoL
    else
        mdot = mdot * rhoR
    end
    
    # p half
    ph = p_p5*pL + p_m5*pR

    return mdot, ph
end

function AUSM_plusup_half(rhoL, rhoR, UL, UR, pL, pR, Minf, g, i,j)
    # param
    beta  = 1/8
    kp    = 0.25
    sigma = 1.0
    ku    = 0.75
    
    # L, R
    aL = (g * pL / rhoL)^0.5
    aR = (g * pR / rhoR)^0.5

    # half
    ah   = 0.5*( aL + aR )
    rhoh = 0.5*( rhoL + rhoR )

    Mbar = (( UL^2 + UR^2 ) / ( 2 * ah^2 ))^0.5
    Mo   = (min(1,max( Mbar^2, Minf^2 )))^0.5
    fa   = Mo * (2-Mo)

    alpha = 3/16 * ( -4 + 5*fa^2 )

    ML = UL/ah
    MR = UR/ah

    M_p4 = 0
    M_m4 = 0
    p_p5 = 0
    p_m5 = 0
    if abs(ML) >= 1
        M_p4 = 0.5*(ML + abs(ML))
        p_p5 = 0.5*(ML + abs(ML)) / ML
    else
        Mtp  = 0.25*(ML + 1)^2
        Mtm  = -0.25*(ML - 1)^2
        M_p4 = Mtp * (1 - 16*beta*Mtm)
        p_p5 = Mtp * ((2-ML) - 16*alpha*ML*Mtm)
    end
    

    if abs(MR) >= 1
        M_m4 = 0.5*(MR - abs(MR))
        p_m5 = 0.5*(MR - abs(MR)) / MR
    else
        Mtp  = 0.25*(MR + 1)^2
        Mtm  = -0.25*(MR - 1)^2
        M_m4 = Mtm * (1 + 16*beta*Mtp)
        p_m5 = Mtm * ((-2-MR) + 16*alpha*MR*Mtp)
    end
    
    # M half
    Mp = -kp * max( 1 - sigma*Mbar^2, 0.0) * (pR - pL)/(rhoh*ah^2)
    
    Mh = M_p4 + M_m4 + Mp/fa
    
    #=
    if i==100 && j==200
        println(" MMM_AUSM+up ")
        println(M_p4)
        println(M_m4)
        println(Mp/fa)
    end
    =#
        
    # mdot half
    mdot = ah * Mh
    if Mh > 0
        mdot = mdot * rhoL
    else
        mdot = mdot * rhoR
    end
    
    # p half
    pu = -ku * p_p5 * p_m5 * (rhoL + rhoR) * ah *(UR - UL)

    ph = p_p5*pL + p_m5*pR + fa * pu
    return mdot, ph
end

function set_Minf(bdcon, specific_heat_ratio, Rd, nval)
    Minf = 0.0
    check = 0
    for i in 1:4
        if Int(bdcon[i][1]) == 0 || Int(bdcon[i][1]) == 5
            rho = bdcon[i][2]
            u = bdcon[i][3]
            v = bdcon[i][4]
            p = bdcon[i][5]

            a = (specific_heat_ratio * p / rho)^0.5
            Minf = (u^2+v^2)^0.5 / a

            println(" Minf = " * string(Minf))

            check = 1

        elseif Int(bdcon[i][1]) == 6
            rho = bdcon[i][2]
            u = bdcon[i][3]
            v = bdcon[i][4]
            T = bdcon[i][nval+2]

            p = (rho*Rd) * T

            a = (specific_heat_ratio * p / rho)^0.5
            Minf = (u^2+v^2)^0.5 / a
            
            println(" Minf = " * string(Minf))

            check = 1
        end
    end

    if check == 0
        println(" Minf error ")
        println("  ")
        println(" if you don't use inlet condition, ")
        println(" don't use AUSM+up ")
        println("  ")
        throw(UndefVarError(:x))
    end

    return Minf
end


function SLAU_half(rhoL, rhoR, UL, UR, pL, pR, velocity, specific_heat_ratio)
    mdot =1.0
    ph =1.0

    # L, R
    aL = (specific_heat_ratio * pL / rhoL)^0.5
    aR = (specific_heat_ratio * pR / rhoR)^0.5

    # half
    ah   = 0.5*( aL + aR )

    ML = UL/ah
    MR = UR/ah

    temp1 = -max( min( ML, 0.0 ), -1.0)
    temp2 =  min( max( MR, 0.0 ), 1.0)
    g = temp1 * temp2

    barUL = (1-g) * (rhoL*abs(UL) + rhoR*abs(UR))/(rhoL + rhoR) + g*abs(UL)
    barUR = (1-g) * (rhoL*abs(UL) + rhoR*abs(UR))/(rhoL + rhoR) + g*abs(UR)

    hatM = min(1.0, 1/ah * velocity)
    chi  = (1-hatM)^2
        
    beta_a = 0
    beta_b = 0
    if abs(MR) >= 1
        beta_a = 0.5 * (1+sign(-MR))
    else
        beta_a = 0.25 * (2+MR) * (MR-1)^2
    end
    if abs(ML) >= 1
        beta_b = 0.5 * (1+sign(ML))
    else
        beta_b = 0.25 * (2-ML) * (ML+1)^2
    end
    
    # mdot half
    mdot = 0.5 * (rhoL*(UL+abs(barUL)) + rhoR*(UR-abs(barUR)) - chi/ah*(pR-pL))
    
    # p half
    ph = 0.5*(pR+pL) + 0.5*(beta_b-beta_a)*(pL-pR) + (1-chi) * (beta_b+beta_a-1) * 0.5*(pR+pL)

    # SLAU2
    #ph = 0.5*(pR+pL) + 0.5*(beta_b-beta_a) + ((UR^2+UL^2)/2)^0.5 * (beta_b+beta_a-1) * 0.5*(rhoR+rhoL)*ah

    return mdot, ph
end
