# src path
src_path = "C:\\Users\\koike\\Desktop\\git\\general_2d_NS_LES\\src\\"


# main (変更しないこと)
src_read="read_grid.jl"
include(src_path*src_read)
src_read="read_para.jl"
include(src_path*src_read)
src_read="output.jl"
include(src_path*src_read)

src_read="advection_term.jl"
include(src_path*src_read)
src_read="boundary.jl"
include(src_path*src_read)
src_read="converge.jl"
include(src_path*src_read)
src_read="misc.jl"
include(src_path*src_read)
src_read="rhs.jl"
include(src_path*src_read)
src_read="value_setup.jl"
include(src_path*src_read)
src_read="turb_val.jl"
include(src_path*src_read)
src_read="viscos_pturb.jl"
include(src_path*src_read)

src_read="allocation.jl"
include(src_path*src_read)
src_read="cal_time_step.jl"
include(src_path*src_read)
src_main="main.jl"
include(src_path*src_main)


