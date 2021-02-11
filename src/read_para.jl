import JSON

# ------------------------------------
# Commenting out the json file
# ------------------------------------
function make_json(PARAMDAT)
    fff = []
    open(PARAMDAT, "r") do f
        fff = read(f,String)
    end

    fff = replace(fff,r"#(.+)\n" => "\n")
    
    read_PARAMDAT = "read_"*PARAMDAT
    open(read_PARAMDAT, "w") do f
        write(f,fff)
    end
    return read_PARAMDAT
end 

# ------------------------------------
# prepare reading the json file
# ------------------------------------
function read_json(read_PARAMDAT)
    dict=1
    open(read_PARAMDAT, "r") do f
        dicttxt = read(f,String)  # file information to string
        dict = JSON.parse(dicttxt)  # parse and transform data
    end
    return dict
end

# ------------------------------------
# read PARAMDAT
# ------------------------------------
function read_para(dict)
    # file name
    out_file_front = dict["out_file_front"]
    out_ext        = dict["out_ext"]

    # restart
    restartnum   = dict["restartnum"]
    restart_file = dict["Restart"] * restartnum * dict["in_ext"]
    restartnum   = Int(parse(Float64,restartnum))
    
    # output
    init_small　=　parse(Float64,dict["init_small"])
    norm_ok   　=　parse(Float64,dict["norm_ok"])
    
    # time step
    time_integ   = dict["time_integ"]
    nt           = Int(parse(Float64,dict["nt"]))
    dt           = parse(Float64,dict["dt"])
    every_outnum = Int(parse(Float64,dict["every_outnum"]))
    in_nt        = Int(parse(Float64,dict["inner_step"]))
    dtau         = parse(Float64,dict["dtau"])
    cfl          = parse(Float64,dict["CFL"])

    # scheme
    ad_scheme = Int(parse(Float64,dict["advection"]))
    
    # Initial values
    init_rho = parse(Float64,dict["init_rho"])
    init_u   = parse(Float64,dict["init_u"])
    init_v   = parse(Float64,dict["init_v"])
    init_p   = parse(Float64,dict["init_p"])
    init_T   = parse(Float64,dict["init_T"])

    # constant    
    specific_heat_ratio = parse(Float64,dict["specific_heat_ratio"])
    Rd = parse(Float64,dict["Rd"])

    # boundary conditions
    #= 
    bdcon[i][j]
    i:bd number (1,2,3..)
    j=1:bd number (0~4)
    j=2:rho
    j=3:u
    j=4:v
    j=5:p
    =#
    bdcon = []
    k = 1
    while true
        try
            bd=[parse(Int,dict["bd"*string(k)*"_con"]),   parse(Float64,dict["bd"*string(k)*"_rho"]),
                parse(Float64,dict["bd"*string(k)*"_u"]), parse(Float64,dict["bd"*string(k)*"_v"]),
                parse(Float64,dict["bd"*string(k)*"_p"]), parse(Float64,dict["bd"*string(k)*"_T"])]
            push!(bdcon,bd)
        catch
            break
        end
        k+=1
    end
    
    return out_file_front, out_ext, restartnum, restart_file, init_small, norm_ok,
            time_integ, nt, dt, every_outnum, in_nt, dtau, cfl, ad_scheme,
            init_rho, init_u, init_v, init_p, init_T, specific_heat_ratio, Rd, bdcon
end

function input_para(PARAMDAT)
    read_PARAMDAT = make_json(PARAMDAT)
    dict          = read_json(read_PARAMDAT)
    out_file_front, out_ext, restartnum, restart_file, init_small, norm_ok,
    time_integ, nt, dt, every_outnum, in_nt, dtau, cfl, ad_scheme,
    init_rho, init_u, init_v, init_p, init_T, specific_heat_ratio, Rd, bdcon = read_para(dict)
    println("fin reading para")

    return out_file_front, out_ext, restartnum, restart_file, init_small, norm_ok,
            time_integ, nt, dt, every_outnum, in_nt, dtau, cfl, ad_scheme,
            init_rho, init_u, init_v, init_p, init_T, specific_heat_ratio, Rd, bdcon
end