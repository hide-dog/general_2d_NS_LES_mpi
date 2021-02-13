function set_boundary(Qbase, cellxmax, cellymax, vecAx, vecAy, bdcon, Rd, g, nval, boundary_number, rank, size, comm)

    divn = zeros(Int, 4)
    temp = Array{Float64}(undef,2)
    
    # 2:imaginary cell
    icell = 2
    send_arrayx = zeros(cellxmax*icell*nval)
    send_arrayy = zeros(cellymax*icell*nval)
    rec_arrayx = zeros(cellxmax*icell*nval)
    rec_arrayy = zeros(cellymax*icell*nval)

    #println(boundary_number)

    # x-
    if boundary_number[1] == -41
        # parallel with communication
        
        ite = 1
        for l in 1:nval
            for j in 1:cellymax
                for i in icell+1:icell+icell  # 周期境界条件
                    send_arrayy[ite] = Qbase[i,j,l]
                    ite += 1
                end
            end
        end

        # 送る
        send_rank = rank + (size^0.5-1) * (size^0.5)
        MPI.Isend(send_arrayy, send_rank, rank+200, comm)
        
        # cellから受け取る
        rec = MPI.Irecv!(rec_arrayy, send_rank, send_rank+100, comm)
        
        # wait 
        MPI.Waitall!([rec])
        
        # 代入
        ite = 1
        for l in 1:nval
            for j in 1:cellymax
                for i in icell+1:icell+icell
                    Qbase[i,j,l] = rec_arrayy[ite]
                    ite += 1
                end
            end
        end

    elseif boundary_number[1] < 0
        Qbase = boundary_condition_xm(Qbase, cellxmax, cellymax, vecAx, vecAy, bdcon, Rd, g, nval)
    else
        # parallel with communication
        
        ite = 1
        for l in 1:nval
            for j in 1:cellymax
                for i in 1:icell
                    send_arrayy[ite] = Qbase[i,j,l]
                    ite += 1
                end
            end
        end

        # south に送る
        MPI.Isend(send_arrayy, boundary_number[2], rank+200, comm)
        
        # north of south cell から受け取る
        rec = MPI.Irecv!(rec_arrayy, boundary_number[1], boundary_number[1]+100, comm)
        
        # wait 
        MPI.Waitall!([rec])
        
        # 代入
        ite = 1
        for l in 1:nval
            for j in 1:cellymax
                for i in 1:icell
                    Qbase[i,j,l] = rec_arrayy[ite]
                    ite += 1
                end
            end
        end

        # use ?
        divn[1] = 1
    end

    # x+
    if  boundary_number[2] == -41
        # parallel with communication
        
        ite = 1
        for l in 1:nval
            for j in 1:cellymax
                for i in cellxmax-icell-(icell+1):cellxmax-icell  # 周期境界条件
                    send_arrayy[ite] = Qbase[i,j,l]
                    ite += 1
                end
            end
        end

        # 送る
        send_rank = rank - (size^0.5-1) * (size^0.5)
        MPI.Isend(send_arrayy, send_rank, rank+200, comm)
        
        # cellから受け取る
        rec = MPI.Irecv!(rec_arrayy, send_rank, send_rank+100, comm)
        
        # wait 
        MPI.Waitall!([rec])
        
        # 代入
        ite = 1
        for l in 1:nval
            for j in 1:cellymax
                for i in cellxmax-icell-(icell+1):cellxmax-icell  # 周期境界条件
                    Qbase[i,j,l] = rec_arrayy[ite]
                    ite += 1
                end
            end
        end

    elseif boundary_number[2] < 0
        Qbase = boundary_condition_xp(Qbase, cellxmax, cellymax, vecAx, vecAy, bdcon, Rd, g, nval)
    else
        # parallel with communication

        ite = 1
        for l in 1:nval
            for j in 1:cellymax
                for i in cellxmax-icell+1:cellxmax
                    send_arrayy[ite] = Qbase[i,j,l]
                    ite += 1
                end
            end
        end

        # south に送る
        MPI.Isend(send_arrayy, boundary_number[1], rank+100, comm)
        
        # north of south cell から受け取る
        rec = MPI.Irecv!(rec_arrayy, boundary_number[2], boundary_number[2]+200, comm)
        
        # wait 
        MPI.Waitall!([rec])
        
        # 代入
        ite = 1
        for l in 1:nval
            for j in 1:cellymax
                for i in cellxmax-icell+1:cellxmax
                    Qbase[i,j,l] = rec_arrayy[ite]
                    ite += 1
                end
            end
        end

        # use ?
        divn[2] = 1
    end

    # y-
    if  boundary_number[3] == -41
        # parallel with communication
        
        ite = 1
        for l in 1:nval
            for j in icell+1:icell+icell  # 周期境界条件
                for i in 1:cellxmax
                    send_arrayy[ite] = Qbase[i,j,l]
                    ite += 1
                end
            end
        end

        # 送る
        send_rank = rank + (size^0.5-1)
        MPI.Isend(send_arrayy, send_rank, rank+200, comm)
        
        # cellから受け取る
        rec = MPI.Irecv!(rec_arrayy, send_rank, send_rank+100, comm)
        
        # wait 
        MPI.Waitall!([rec])
        
        # 代入
        ite = 1
        for l in 1:nval
            for j in icell+1:icell+icell  # 周期境界条件
                for i in 1:cellxmax
                    Qbase[i,j,l] = rec_arrayy[ite]
                    ite += 1
                end
            end
        end

    elseif boundary_number[3] < 0
        Qbase = boundary_condition_ym(Qbase, cellxmax, cellymax, vecAx, vecAy, bdcon, Rd, g, nval)
    else
        # parallel with communication
        
        ite = 1
        for l in 1:nval
            for j in 1:icell
                for i in 1:cellymax
                    send_arrayy[ite] = Qbase[i,j,l]
                    ite += 1
                end
            end
        end

        # south に送る
        MPI.Isend(send_arrayy, boundary_number[4], rank+400, comm)
        
        # north of south cell から受け取る
        rec = MPI.Irecv!(rec_arrayy, boundary_number[3], boundary_number[3]+300, comm)
        
        # wait 
        MPI.Waitall!([rec])
        
        ite = 1
        for l in 1:nval
            for j in 1:icell
                for i in 1:cellxmax
                    Qbase[i,j,l] = send_arrayy[ite]
                    ite += 1
                end
            end
        end

        # use ?
        divn[3] = 1
    end

    # y+
    if boundary_number[4] == -41
        
        ite = 1
        for l in 1:nval
            for j in cellymax-icell-(icell+1):cellymax-icell  # 周期境界条件
                for i in 1:cellxmax
                    send_arrayy[ite] = Qbase[i,j,l]
                    ite += 1
                end
            end
        end

        # 送る
        send_rank = rank - (size^0.5-1)
        MPI.Isend(send_arrayy, send_rank, rank+200, comm)
        
        # cellから受け取る
        rec = MPI.Irecv!(rec_arrayy, send_rank, send_rank+100, comm)
        
        # wait 
        MPI.Waitall!([rec])
        
        # 代入
        ite = 1
        for l in 1:nval
            for j in cellymax-icell-(icell+1):cellymax-icell  # 周期境界条件
                for i in 1:cellxmax
                    Qbase[i,j,l] = rec_arrayy[ite]
                    ite += 1
                end
            end
        end
    elseif boundary_number[4] < 0
        Qbase = boundary_condition_yp(Qbase, cellxmax, cellymax, vecAx, vecAy, bdcon, Rd, g, nval)
    else
        # parallel with communication
        
        ite = 1
        for l in 1:nval
            for j in cellymax-icell+1:cellymax
                for i in 1:cellymax
                    send_arrayy[ite] = Qbase[i,j,l]
                    ite += 1
                end
            end
        end

        # south に送る
        MPI.Isend(send_arrayy, boundary_number[3], rank+300, comm)
        
        # north of south cell から受け取る
        rec = MPI.Irecv!(rec_arrayy, boundary_number[4], boundary_number[4]+400, comm)
        
        # wait 
        MPI.Waitall!([rec])
        
        ite = 1
        for l in 1:nval
            for j in cellymax-icell+1:cellymax
                for i in 1:cellxmax
                    Qbase[i,j,l] = send_arrayy[ite]
                    ite += 1
                end
            end
        end

        # use ?
        divn[4] = 1
    end

    return Qbase
end

# ------------------------------------
# boundary conditions x- (x[1] & x[2])
# ------------------------------------
function boundary_condition_xm(Qbase, cellxmax, cellymax, vecAx, vecAy, bdcon, Rd, g, nval)
    """
    bdcon[i][j]
    i : bd number(x-, x+, y-, y+)
    j=1-6 : "bd1_con":"2",
            "bd1_rho":"1.0",
            "bd1_u"  :"300.0",
            "bd1_v"  :"0.0",
            "bd1_p"  :"1.0",
            "bd1_T"  :"300.0",
    """

    # bd1 = x-
    if Int(bdcon[1][1]) == 11
        for j in 1:cellymax
            for l in 1:nval
                Qbase[1,j,l] = bdcon[1][l+1]
            end
        end
    elseif Int(bdcon[1][1]) == 21
        for j in 1:cellymax
            for l in 1:nval
                Qbase[1,j,l] = Qbase[2,j,l]
            end
        end
    elseif Int(bdcon[1][1]) == 31
        for j in 1:cellymax
            for l in 1:nval
                Qbase[1,j,l] = Qbase[2,j,l]
            end

            xvec = vecAx[2,j,1]
            yvec = vecAx[2,j,2]
            u = Qbase[2,j,2]
            v = Qbase[2,j,3]

            Qbase[1,j,2] = ((-xvec^2+yvec^2)*u-2*xvec*yvec*v)/(xvec^2+yvec^2)
            Qbase[1,j,3] = (-2*xvec*yvec*u+(xvec^2-yvec^2)*v)/(xvec^2+yvec^2)            
        end
    elseif Int(bdcon[1][1]) == 32
        for j in 1:cellymax
            for l in 1:nval
                Qbase[1,j,l] = Qbase[2,j,l]
            end

            u = Qbase[2,j,2]
            v = Qbase[2,j,3]

            Qbase[1,j,2] = -u
            Qbase[1,j,3] = -v
        end
    elseif Int(bdcon[1][1]) == 41
        for j in 1:cellymax
            for l in 1:nval
                Qbase[1,j,l] = Qbase[cellxmax-1,j,l]
            end
        end
    else
        println("------------------------")
        println(" boundary error ")
        println(" please check boundary number ")
        println("------------------------")
        throw(UndefVarError(:x))
    end

    return Qbase
end

# ------------------------------------
# boundary conditions x+ (x[cellxmax-1] & x[cellxmax])
# ------------------------------------
function boundary_condition_xp(Qbase, cellxmax, cellymax, vecAx, vecAy, bdcon, Rd, g, nval)
    
    # bd2 = x+
    if Int(bdcon[2][1]) == 11
        for j in 1:cellymax
            for l in 1:nval
                Qbase[cellxmax,j,l] = bdcon[2][l+1]
            end
        end
    elseif Int(bdcon[2][1]) == 21
        for j in 1:cellymax
            for l in 1:nval
                Qbase[cellxmax,j,l] = Qbase[cellxmax-1,j,l]
            end
        end
    elseif Int(bdcon[2][1]) == 31
        for j in 1:cellymax
            for l in 1:nval
                Qbase[cellxmax,j,l] = Qbase[cellxmax-1,j,l]
            end

            xvec = vecAx[cellymax,j,1]
            yvec = vecAx[cellymax,j,2]
            u = Qbase[cellxmax-1,j,2]
            v = Qbase[cellxmax-1,j,3]

            Qbase[cellxmax,j,2] = ((-xvec^2+yvec^2)*u-2*xvec*yvec*v)/(xvec^2+yvec^2)
            Qbase[cellxmax,j,3] = (-2*xvec*yvec*u+(xvec^2-yvec^2)*v)/(xvec^2+yvec^2)
        end
    elseif Int(bdcon[2][1]) == 32
        for j in 1:cellymax
            for l in 1:nval
                Qbase[cellxmax,j,l] = Qbase[cellxmax-1,j,l]
            end
            u = Qbase[cellxmax-1,j,2]
            v = Qbase[cellxmax-1,j,3]

            Qbase[cellxmax,j,2] = -u
            Qbase[cellxmax,j,3] = -v
        end
    elseif Int(bdcon[2][1]) == 41
        for j in 1:cellymax
            for l in 1:nval
                Qbase[cellxmax,j,l] = Qbase[2,j,l]
            end
        end
    else
        println("------------------------")
        println(" boundary error ")
        println(" please check boundary number ")
        println("------------------------")
        throw(UndefVarError(:x))
    end

    return Qbase
end

# ------------------------------------
# boundary conditions y- (y[1] & y[2])
# ------------------------------------
function boundary_condition_ym(Qbase, cellxmax, cellymax, vecAx, vecAy, bdcon, Rd, g, nval)
    
    # bd3 = y-
    if Int(bdcon[3][1]) == 11
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,1,l] = bdcon[3][l+1]
            end
        end
    elseif Int(bdcon[3][1]) == 21
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,1,l] = Qbase[i,2,l]
            end
        end
    elseif Int(bdcon[3][1]) == 31
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,1,l] = Qbase[i,2,l]
            end

            xvec = vecAy[i,2,1]
            yvec = vecAy[i,2,2]
            u = Qbase[i,2,2]
            v = Qbase[i,2,3]
            
            Qbase[i,1,2] = ((-xvec^2+yvec^2)*u-2*xvec*yvec*v)/(xvec^2+yvec^2)
            Qbase[i,1,3] = (-2*xvec*yvec*u+(xvec^2-yvec^2)*v)/(xvec^2+yvec^2)
        end
    elseif Int(bdcon[3][1]) == 32
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,1,l] = Qbase[i,2,l]
            end

            u = Qbase[i,2,2]
            v = Qbase[i,2,3]

            Qbase[i,1,2] = -u
            Qbase[i,1,3] = -v
        end
    elseif Int(bdcon[3][1]) == 41
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,1,l] = Qbase[i,cellymax-1,l]
            end 
        end
    elseif Int(bdcon[3][1]) == 33
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,1,l] = Qbase[i,2,l]
            end

            u = Qbase[i,2,2]
            v = Qbase[i,2,3]

            Qbase[i,1,2] = -u
            Qbase[i,1,3] = -v

            Tw = bdcon[3][nval+2]
            p  = Qbase[i,1,4]
            Qbase[i,1,1] = p/(Tw*Rd)
        end
    else
        println("------------------------")
        println(" boundary error ")
        println(" please check boundary number ")
        println("------------------------")
        throw(UndefVarError(:x))
    end

    return Qbase
end

# ------------------------------------
# boundary conditions y+ (y[cellymax-1] & y[cellymax])
# ------------------------------------
function boundary_condition_yp(Qbase, cellxmax, cellymax, vecAx, vecAy, bdcon, Rd, g, nval)
    
    # bd4 = y+
    if Int(bdcon[4][1]) == 11
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,cellymax,l] = bdcon[4][l+1]
            end
        end
    elseif Int(bdcon[4][1]) == 21
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,cellymax,l] = Qbase[i,cellymax-1,l]
            end
        end
    elseif Int(bdcon[4][1]) == 31
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,cellymax,l] = Qbase[i,cellymax-1,l]
            end

            xvec = vecAy[i,cellymax,1]
            yvec = vecAy[i,cellymax,2]
            u = Qbase[i,cellymax-1,2]
            v = Qbase[i,cellymax-1,3]

            Qbase[i,cellymax,2] = ((-xvec^2+yvec^2)*u-2*xvec*yvec*v)/(xvec^2+yvec^2)
            Qbase[i,cellymax,3] = (-2*xvec*yvec*u+(xvec^2-yvec^2)*v)/(xvec^2+yvec^2)
        end
    elseif Int(bdcon[4][1]) == 32
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,cellymax,l] = Qbase[i,cellymax-1,l]
            end

            u = Qbase[i,cellymax-1,2]
            v = Qbase[i,cellymax-1,3]

            Qbase[i,cellymax,2] = -u
            Qbase[i,cellymax,3] = -v
        end
    elseif Int(bdcon[4][1]) == 41
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,cellymax,l] = Qbase[i,2,l]
            end 
        end
    elseif Int(bdcon[4][1]) == 101
        temp1 = Int64(round(cellxmax/4)+1) # 1/4流出
        temp2 = Int64((temp1-1)*3+1)       # 1/4~3/4流入
        for i in 1:temp1
            for l in 1:nval
                Qbase[i,cellymax,l] = Qbase[i,cellymax-1,l]
            end
        end
        for i in temp1+1:temp2
            for l in 1:nval
                Qbase[i,cellymax,l] = bdcon[4][l+1]
            end
        end
        for i in temp2+1:cellxmax
            for l in 1:nval
                Qbase[i,cellymax,l] = Qbase[i,cellymax-1,l]
            end
        end
    elseif Int(bdcon[4][1]) == 102
        temp1 = Int64(round(cellxmax/4)+1) # 1/4流出
        temp2 = Int64((temp1-1)*3+1)       # 1/4~3/4流入
        for i in 1:temp1
            for l in 1:nval
                Qbase[i,cellymax,l] = Qbase[i,cellymax-1,l]
            end
        end
        for i in temp1+1:temp2
            for l in 1:nval
                Qbase[i,cellymax,l] = bdcon[4][l+1]
            end
            T = bdcon[4][nval+2]
            p = (Qbase[i,cellymax,1]*Rd) * T
            Qbase[i,cellymax,4] = p
        end
        for i in temp2+1:cellxmax
            for l in 1:nval
                Qbase[i,cellymax,l] = Qbase[i,cellymax-1,l]
            end
        end
    elseif Int(bdcon[4][1]) == 12
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,cellymax,l] = bdcon[4][l+1]
            end
            T = bdcon[4][nval+2]
            p = (Qbase[i,cellymax,1]*Rd) * T
            Qbase[i,cellymax,4] = p
        end
    elseif Int(bdcon[4][1]) == 9
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,cellymax,l] = bdcon[4][l+1]
            end

            Ayx = vecAy[i,cellymax,1] / (vecAy[i,cellymax,1]^2 + vecAy[i,cellymax,2]^2)^0.5
            Ayy = vecAy[i,cellymax,2] / (vecAy[i,cellymax,1]^2 + vecAy[i,cellymax,2]^2)^0.5
            
            #=
            Un = Qbase[i,cellymax-1,2]*Ayx + Qbase[i,cellymax-1,2]*Ayy
            if Un < 0

            end
            =#
        end
    else
        println("------------------------")
        println(" boundary error ")
        println(" please check boundary number ")
        println("------------------------")
        throw(UndefVarError(:x))
    end

    return Qbase
end

# ------------------------------------
# check boundary condition
# ------------------------------------
function check_bd(bdcon)
    bdnum = [11, 12, 101, 102, 21, 31, 32, 33, 41]
    for l in 1:4
        di = 0
        for nn in bdnum
            if Int(bdcon[l][1]) == nn
                di = 1 
            end
        end

        if di == 0
            println("\n check boundary condition ! \n")
            throw(UndefVarError(:x))
        end
    end
end

# ------------------------------------
# set boundary condition number for mpi
# ------------------------------------
function set_boundary_number(bdcon, rank, size)
    #=
    bd_NSEW = [N, S, E, W]
    北南東西に対応する境界条件を入れる．
    number <  0 : 境界条件,bdcon
    number >= 0 : 対応するセルのrank
    =#
    #=
    evaluation system

    parallel num = 4

      rank      div(rank, size^0.5)    rem(rank, size^0.5)

      2 | 3            1 | 1                  0 | 1
     -------          -------                -------
      0 | 1            0 | 0                  0 | 1
    
    parallel num = 16

          rank              div(rank, size^0.5)       rem(rank, size^0.5)

    12 | 13 | 14 | 15          3  | 3 | 3  | 3          0  | 1 | 2  | 3  
    ------------------        ------------------       ------------------
    8  | 9  | 10 | 11          2  | 2 | 2  | 2          0  | 1 | 2  | 3  
    ------------------        ------------------       ------------------
    4  | 5  | 6  | 7           1  | 1 | 1  | 1          0  | 1 | 2  | 3  
    ------------------        ------------------       ------------------
    0  | 1  | 2  | 3           0  | 0 | 0  | 0          0  | 1 | 2  | 3  

    0 : upper or left edge
    sqrt(parallel)-1 : lower or right edge 

    =#
    bd_NSEW = zeros(Int, 4)
    y = div(rank, size^0.5)
    x = rem(rank, size^0.5)

    if y == 0 && y == Int(size^0.5 - 1)
        bd_NSEW[1] = - Int(bdcon[2][1])
        bd_NSEW[2] = - Int(bdcon[1][1])
    elseif y == 0
        bd_NSEW[1] = rank + size
        bd_NSEW[2] = - Int(bdcon[1][1])
    elseif y == Int(size^0.5 - 1)
        bd_NSEW[1] = - Int(bdcon[2][1])
        bd_NSEW[2] = rank - size
    else
        bd_NSEW[1] = rank + size
        bd_NSEW[2] = rank - size
    end

    if x == 0 && x == Int(size^0.5 - 1)
        bd_NSEW[3] = - Int(bdcon[4][1])
        bd_NSEW[4] = - Int(bdcon[3][1])
    elseif x == 0
        bd_NSEW[3] = rank + 1
        bd_NSEW[4] = - Int(bdcon[3][1])
    elseif x == Int(size^0.5 - 1)
        bd_NSEW[3] = - Int(bdcon[4][1])
        bd_NSEW[4] = rank - 1
    else
        bd_NSEW[3] = rank + 1
        bd_NSEW[4] = rank - 1
    end
    return bd_NSEW
end