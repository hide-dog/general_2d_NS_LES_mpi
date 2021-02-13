import MPI

# ------------------------------------
# mpi
# ------------------------------------
function set_mpi(cellxmax, cellymax)
    comm = MPI.COMM_WORLD        # MPIプロセスの呼び出し
    rank = MPI.Comm_rank(comm)   # MPIの分割, 各プロセス番号(0~n)
    size = MPI.Comm_size(comm)   # MPIの分割, 各プロセス番号(0~n)
    
    num_div = 0.0      # 各方向の分割数
    if floor(size^0.5 * 10000) % 10000 != 0 || size <= 0.0
        println("------------------------\n")
        println(" mpi error \n")
        println(" please change parallel number to be an number of squares \n")
        println("------------------------")
        throw(UndefVarError(:x))
    else
        num_div = Int(size^0.5)
    end

    Nx = cellxmax-4 # 仮想セル分を引く
    Ny = cellymax-4 # 仮想セル分を引く

    #=
    並列時のstart i, end i, 分割数cellNxを定義
    =#
    cellNx, ista, iend, cellNy, jsta, jend = start_and_end0(Nx, Ny, num_div, rank)

    return cellNx, ista, iend, cellNy, jsta, jend, num_div, rank, size, comm
end

# ---------------------------------------------
# mpi division
# ---------------------------------------------
function start_and_end0(Nx, Ny, num_div, rank)

    if Nx % num_div != 0
        println("------------------------\n")
        println(" mpi error \n ")
        println(" please change Nx \n ")
        println(" Nx = ")
        println(Nx)
        println(" \n ")
        println(" num_div = ")
        println(num_div)
        println("------------------------")
    elseif Ny % num_div != 0
        println("------------------------\n")
        println(" mpi error \n ")
        println(" please change Ny \n ")
        println(" Ny = ")
        println(Ny)
        println(" \n ")
        println(" num_div = ")
        println(num_div)
        println("------------------------")
    end

    cellNx = div(Nx, num_div)   # プロセスの分割数
    cellNy = div(Ny, num_div)   # プロセスの分割数
    ista = rank * cellNx + 1  # プロセス毎のstart i (0-> 0*n+1 = 1,        1-> 1*n+1 = n+1)
    iend = ista + cellNx - 1  # プロセス毎のend i   (0-> 0*n+1 + n-1 = n , 1-> 1*n+1 + n-1 = 2n)
    jsta = rank * cellNy + 1  # プロセス毎のstart i (0-> 0*n+1 = 1,        1-> 1*n+1 = n+1)
    jend = ista + cellNy - 1  # プロセス毎のend i   (0-> 0*n+1 + n-1 = n , 1-> 1*n+1 + n-1 = 2n)
    
    return cellNx, ista, iend, cellNy, jsta, jend
end

function set_value_for_mpi(Qbase, cellx_mpi, celly_mpi, volume, cellcenter, wally, dx, dy,
                        Qbase_all, volume_all, cellcenter_all, wally_all, dx_all, dy_all, nval, num_div, rank)
    
    is = 3
    ie = cellx_mpi-2
    js = 3
    je = celly_mpi-2
    for l in 1:nval
        for i in is:ie
            for j in js:je
                ii, jj = Conversion_to_coordinate_numbers(i, j, cellx_mpi, celly_mpi, num_div, rank)
                Qbase[i,j,l]  = Qbase_all[ii, jj, l]
            end
        end
    end
    for i in is:ie
        for j in js:je
            ii, jj = Conversion_to_coordinate_numbers(i, j, cellx_mpi, celly_mpi, num_div, rank)
            volume[i,j] = volume_all[ii, jj]
            wally[i,j]  = wally_all[ii, jj]
            dx[i,j]     = dx_all[ii, jj]
            dy[i,j]     = dy_all[ii, jj]
        end
    end
    
    return Qbase, volume, cellcenter, wally, dx, dy
end

function Conversion_to_coordinate_numbers(i, j, cellx_mpi, celly_mpi, num_div, rank)
    #=
    div system

      rank      div(rank, size^0.5)    rem(rank, size^0.5)

      2 | 3            1 | 1                  0 | 1
     -------          -------                -------
      0 | 1            0 | 0                  0 | 1
      
      y
      ↑
       →x　

    =#
    rankx = rem(rank, num_div)
    ranky = div(rank, num_div)

    ii = i + rankx*(cellx_mpi-2)
    jj = j + ranky*(celly_mpi-2)

    return ii, jj    
end

function set_vecA_mpi(vecAx_all, vecAy_all, vecAx, vecAy, cellx_mpi, celly_mpi, num_div, rank)
    #=
    vecA system

                    vecAy[2,3]
                |---------------|
                |               |
                |               |
    vecAx[2,2]  |   cell(2,2)   | vecAx[3,2]
                |               |
                |               |
                |---------------|
                    vecAy[2,2]

    =#
    is = 1
    ie = cellx_mpi+1
    js = 1
    je = celly_mpi
    for l in 1:2
        for i in is:ie
            for j in js:je
                ii, jj = Conversion_to_coordinate_numbers_atbound(i, j, cellx_mpi, celly_mpi, num_div, rank)
                vecAx[i,j,l]  = vecAx_all[ii, jj, l]
            end
        end
    end

    is = 1
    ie = cellx_mpi
    js = 1
    je = celly_mpi+1
    for l in 1:2
        for i in is:ie
            for j in js:je
                ii, jj = Conversion_to_coordinate_numbers_atbound(i, j, cellx_mpi, celly_mpi, num_div, rank)
                vecAy[i,j,l]  = vecAy_all[ii, jj, l]
            end
        end
    end

    return vecAx, vecAy
end

function Conversion_to_coordinate_numbers_atbound(i, j, cellx_mpi, celly_mpi, num_div, rank)
    #=
    div system

      rank      div(rank, size^0.5)    rem(rank, size^0.5)

      2 | 3            1 | 1                  0 | 1
     -------          -------                -------
      0 | 1            0 | 0                  0 | 1

    =#
    rankx = rem(rank, num_div)
    ranky = div(rank, num_div)

    ii = i + rankx*(cellx_mpi)
    jj = j + ranky*(celly_mpi)

    return ii, jj    
end

function mpi_div_3d(ite, ls, le, js, je, is, ie)
    #= ite: iteration number (ista, iend generated by start_and_end(N,size))
        ls: start number at l 
        le: end number at l 
        js:
        je:
        is:
        ie:
        div():商
        rem():余り
    =#
    Ni = ie - is + 1 # number of i ite
    Nj = je - js + 1 # number of j ite
    #Nl = le - ls + 1 # number of l ite

    l = ls + div(ite-1, Ni*Nj)
    
    ite = ite - div(ite-1, Ni*Nj) * Ni*Nj

    j = js + div(ite-1, Ni)
    
    ite = ite - div(ite-1, Ni) * Ni

    i = is + div(ite-1, 1)

    return i, j, l
end

function mpi_div_2d(ite, js, je, is, ie)
    #= ite: iteration number (ista, iend generated by start_and_end(N,size))
        js:
        je:
        is:
        ie:
        div():商
        rem():余り
    =#
    Ni = ie - is + 1 # number of i ite

    j = js + div(ite-1, Ni)
    
    ite = ite - div(ite-1, Ni) * Ni

    i = is + div(ite-1, 1)

    return i, j
end
