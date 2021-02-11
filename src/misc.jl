# ------------------------------------
# set initial conditions
# ------------------------------------
function set_initQbase(Qbase, cellxmax, cellymax, restart_file, init_rho, init_u, init_v, init_p, init_T,
                    specific_heat_ratio, out_file_front, out_ext, out_dir, restartnum, Rd, nval)

    restart_check = 0
    try Qbase = setup_restart_value(Qbase, cellxmax, cellymax, out_dir, restart_file, nval)
        println("Restart "*restart_file)
        restart_check = 2
    catch 
        restart_check = 1
    end

    if restart_check == 1
        Qbase = setup_init_value(Qbase, cellxmax, cellymax, init_rho, init_u, init_v, init_p)
        println("Start Initial condition")
        restartnum = 0
        output_result(0, Qbase, cellxmax, cellymax, specific_heat_ratio, out_file_front, out_ext, out_dir, Rd, nval)
    end

    return Qbase, restartnum
end

function setup_init_value(Qbase, cellxmax, cellymax, init_rho, init_u, init_v, init_p)
    for j in 1:cellymax
        for i in 1:cellxmax
            Qbase[i,j,1] = init_rho
            Qbase[i,j,2] = init_u
            Qbase[i,j,3] = init_v
            Qbase[i,j,4] = init_p
        end
    end
    return Qbase
end

function setup_restart_value(Qbase, cellxmax, cellymax, out_dir, restart_file, nval)
    skipnum = 1
    fff = []
    open("result/"*restart_file, "r") do f
        fff = read(f, String)
    end 
    fff = split(fff,"\n", keepempty = false)
    
    for i in 1+skipnum:length(fff)
        fff[i] = replace(fff[i]," \r" => "")
    end
    
    k = 1
    for i in 2:cellxmax-1
        for j in 2:cellymax-1
            temp = split(fff[k+skipnum]," ")
            for l in 1:nval
                Qbase[i,j,l] = parse(Float64,temp[l]) 
            end
            k = k+1
        end
    end
    return Qbase
end

# ------------------------------------
# find out if results were diverge
# ------------------------------------
function check_divrege(Qbase, cellxmax, cellymax, Rd, fwrite)
    ite = 0
    for j in 2:cellymax-1
        for i in 2:cellxmax-1        
            T = Qbase[i,j,4]/(Qbase[i,j,1]*Rd)
            if T < 0 || isequal(T, NaN) == true
                open( fwrite, "a" ) do f
                    ai  = @sprintf("%4.0f", i)
                    aj  = @sprintf("%4.0f", j)
                    rho = @sprintf("%4.0e", Qbase[i,j,1])
                    p   = @sprintf("%4.0e", Qbase[i,j,4])

                    write(f, "\n")
                    write(f, " diverge ")
                    write(f, "\n")
                    write(f, " i = "*ai)
                    write(f, "\n")
                    write(f, " j = "*aj)
                    write(f, "\n")
                    write(f, " rho = "*rho)
                    write(f, "\n")
                    write(f, " p = "*p)
                    write(f, "\n")

                end
                ite = 1
            end
        end
    end

    if ite == 1
        println("\n")
        println("\n")
        println(" T<0 ")
        println(" diverge ")
        println("\n")
        throw(UndefVarError(:x))
    end
end

# ------------------------------------
# quicksort
# ------------------------------------
function quicksort(list, first, last, num_cell_data, standard_num)
    # list = zezros(m,n)
    x = list[Int(floor((first+last)/2)),standard_num]
    i = first
    j = last
    
    while true
        while list[i,standard_num] < x
            i = i+1
        end
        while x < list[j,standard_num]
            j = j-1
        end
        
        if (i >= j) 
            break
        else
            for k in 1:num_cell_data
                t = list[i,k]
                list[i,k] = list[j,k]
                list[j,k] = t
            end
            i = i+1
            j = j-1
        end
    end
    
    if (first < i - 1) 
        list = quicksort(list, first, i - 1,num_cell_data,standard_num)
    end
    if (j + 1 < last) 
        list = quicksort(list, j + 1, last,num_cell_data,standard_num)
    end
    return list
end