# ------------------------------------
# calculate right hand side
# ------------------------------------
function setup_RHS(RHS, cellxmax, cellymax, E_adv_hat, F_adv_hat, E_vis_hat, F_vis_hat, nval, volume)
    Threads.@threads for l in 1:nval
        for j in 2:cellymax-1
            for i in 2:cellxmax-1
                RHS[i,j,l] = - (E_adv_hat[i+1,  j, l] - E_vis_hat[i+1,  j, l]) +
                             (E_adv_hat[  i,  j, l] - E_vis_hat[  i,  j, l]) -
                             (F_adv_hat[  i,j+1, l] - F_vis_hat[  i,j+1, l]) +
                             (F_adv_hat[  i,  j, l] - F_vis_hat[  i,  j, l])
            end
        end
    end
    return RHS
end