# ------------------------------------
# calculate Residuals
# ------------------------------------
function set_res(Delta_Qcon_hat, Delta_Qcon_hat_temp, cellxmax, cellymax, nval)
    res = zeros(cellxmax, cellymax, nval)
    for l in 1:nval
        for j in 2:cellymax-1
            for i in 2:cellxmax-1
                res[i,j,l] = abs(Delta_Qcon_hat[i,j,l]-Delta_Qcon_hat_temp[i,j,l])
            end
        end
    end 
    return res
end

# ------------------------------------
# calculate Residuals by norm-2
# ------------------------------------
function check_converge(res, RHS, cellxmax, cellymax, init_small, nval)
    norm2 = zeros(nval)

    tempAxb = zeros(nval)
    tempb = zeros(nval)
    for l in 1:nval
        for j in 2:cellymax-1
            for i in 2:cellxmax-1            
                tempAxb[l] = tempAxb[l] + res[i,j,l]^2
                tempb[l] = tempb[l] + RHS[i,j,l]^2
            end
        end
    end
    for l in 1:nval
        norm2[l] = (tempAxb[l]/(tempb[l]+init_small)) ^ 0.5
    end
    return norm2
end