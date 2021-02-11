# ------------------------------------
# read noumber of nodes
# ------------------------------------
function read_nodenum(skipnum)
    """ 
    xmax: number of nodes in xi direction
    ymax: number of nodes in eta direction
    """
    fff = []
    open("grid/nodesnum", "r") do f
        fff = read(f,String)
    end 
    fff = split(fff,"\n",keepempty=false)
    num_nodes = length(fff) - skipnum

    for i in 1+skipnum:length(fff)
        fff[i] = replace(fff[i]," \r" => "")
    end

    temp = split(fff[2]," ")
    xmax = parse(Int64,temp[1]) 
    ymax = parse(Int64,temp[2]) 
    
    return xmax, ymax
end

# ------------------------------------
# read nodes
# ------------------------------------
function read_nodes(skipnum,xmax,ymax)
    fff = []
    open("grid/nodes", "r") do f
        fff = read(f,String)
    end 
    fff = split(fff,"\n",keepempty=false)
    num_nodes = length(fff) - skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i] = replace(fff[i]," \r" => "")
    end

    nodes = zeros(xmax,ymax,2)
    for i in 1:num_nodes
        temp = split(fff[i+skipnum]," ")

        xnum = parse(Int64,temp[1])
        ynum = parse(Int64,temp[2])
        nodes[xnum,ynum,1] = parse(Float64,temp[3])
        nodes[xnum,ynum,2] = parse(Float64,temp[4]) 
    end
    return nodes
end 

# ------------------------------------
# read vecAx and vecAy
# ------------------------------------
function read_vecA(skipnum,xmax,ymax)
    fff = []
    open("grid/vecAx", "r") do f
        fff = read(f,String)
    end 
    fff = split(fff,"\n",keepempty=false)
    num_cell = length(fff) - skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i] = replace(fff[i]," \r" => "")
    end

    vecAx = zeros(xmax, ymax-1, 2)
    for i in 1:num_cell
        temp = split(fff[i+skipnum]," ")

        x = parse(Int64,temp[1])
        y = parse(Int64,temp[2])
        vecAx[x,y,1] = parse(Float64,temp[3])
        vecAx[x,y,2] = parse(Float64,temp[4]) 
    end

    fff = []
    open("grid/vecAy", "r") do f
        fff = read(f,String)
    end 
    fff = split(fff,"\n",keepempty=false)
    num_cell = length(fff) - skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i] = replace(fff[i]," \r" => "")
    end

    vecAy = zeros(xmax-1, ymax, 2)
    for i in 1:num_cell
        temp = split(fff[i+skipnum]," ")

        x = parse(Int64,temp[1])
        y = parse(Int64,temp[2])
        vecAy[x,y,1] = parse(Float64,temp[3])
        vecAy[x,y,2] = parse(Float64,temp[4]) 
    end
    return vecAx, vecAy
end 

function read_allgrid()
    skip = 1
    xmax, ymax   = read_nodenum(skip)
    nodes        = read_nodes(skip, xmax, ymax)
    vecAx, vecAy = read_vecA(skip, xmax, ymax)
    println("fin reading grid")
    return xmax, ymax, nodes, vecAx, vecAy
end
