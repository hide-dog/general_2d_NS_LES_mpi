function set_wally(nodes, bdcon, wally, cellcenter, cellxmax, cellymax)
	# serch wall
	swith_wall = zeros(4)   # x-, x+, y-, y+を検索する
							# 1:壁
	ite = 0
	for i in 1:4
		if Int(bdcon[i][1]) == 31
			swith_wall[i] = 1
			ite = 1
		elseif Int(bdcon[i][1]) == 32
			swith_wall[i] = 1
			ite = 1
		elseif Int(bdcon[i][1]) == 33
			swith_wall[i] = 1
			ite = 1
		end
	end

	if ite == 0
		println("\n check wall boundary condition ! \n")
		println("\n error at turbulence \n")
		throw(UndefVarError(:x))
	end

	# cal number of ponits
	nop = 0
	if swith_wall[1] == 1
		nop = nop + (cellymax-1)
	end
	if swith_wall[2] == 1
		nop = nop + (cellymax-1)
	end
	if swith_wall[3] == 1
		nop = nop + (cellxmax-1)
	end
	if swith_wall[4] == 1
		nop = nop + (cellxmax-1)
	end

	wallpoint = zeros(nop, 2) # (格子数, 座標)
	ite = 1

	# x-
	if swith_wall[1] == 1
		for j in 2:cellymax-1
			wallpoint[ite,1] = nodes[2,j,1]
			wallpoint[ite,2] = nodes[2,j,2]
			ite = ite + 1
		end
	end
	# x+
	if swith_wall[2] == 1
		for j in 2:cellymax-1
			wallpoint[ite,1] = nodes[2,j,1]
			wallpoint[ite,2] = nodes[2,j,2]
			ite = ite + 1
		end
	end
	# y-
	if swith_wall[3] == 1
		for i in 2:cellxmax-1
			wallpoint[ite,1] = nodes[i,2,1]
			wallpoint[ite,2] = nodes[i,2,2]
			ite = ite + 1
		end
	end
	# y+
	if swith_wall[4] == 1
		for i in 2:cellxmax-1
			wallpoint[ite,1] = nodes[i,2,1]
			wallpoint[ite,2] = nodes[i,2,2]
			ite = ite + 1
		end
	end

	# 距離格納用
	distance_list = zeros(nop,2)
	distance = 0

	# 最大値を2つ抽出
	for j in 2:cellymax-1
		for i in 2:cellxmax-1
			for n in 1:nop
				distance = ((cellcenter[i,j,1]-wallpoint[n,1])^2 +
							(cellcenter[i,j,2]-wallpoint[n,2])^2
							)^0.5
				distance_list[n,1] = n
				distance_list[n,2] = distance
			end

			qdistance_list = copy(distance_list)
			qdistance_list = quicksort(qdistance_list,1,nop,2,2)

			# wall_point
			point1 = Int(qdistance_list[1,1])
			point2 = Int(qdistance_list[2,1])
			
			dis_bottom = ((wallpoint[point1, 1] - wallpoint[point2, 1])^2 + (wallpoint[point1, 2] - wallpoint[point2, 2])^2)^0.5
			
			# O型格子の関係で同じ点が二つあるため
			if dis_bottom == 0.0
				point1 = Int(qdistance_list[1,1])
				point2 = Int(qdistance_list[3,1])
			
				dis_bottom = ((wallpoint[point1, 1] - wallpoint[point2, 1])^2 + (wallpoint[point1, 2] - wallpoint[point2, 2])^2)^0.5
			end

			dis1 = distance_list[point1, 2]
			dis2 = distance_list[point2, 2]
			
			# ヘロンの公式から距離を算出
			s = (dis_bottom + dis1 + dis2)*0.5
			surface = (s*(s-dis_bottom)*(s-dis1)*(s-dis2))^0.5
			wally[i,j] = surface/dis_bottom * 2
		end
	end

	# 仮想セルへの代入
	for i in 2:cellxmax-1
		wally[i,1] = wally[i,2]
		wally[i,cellymax] = wally[i,cellymax-1]
	end
	for j in 2:cellymax-1
		wally[1,j] = wally[2,j]
		wally[cellxmax,j] = wally[cellxmax-1,j]
	end

	return wally, swith_wall
end

function cal_yplus(yplus, Qbase, wally, swith_wall, mu, cellxmax, cellymax, vecAx, vecAy, volume)
	tvAxx = 0.0
    tvAxy = 0.0
    Axx = 0.0
    Axy = 0.0
    tvAyx = 0.0
    tvAyy = 0.0
    Ayx = 0.0
    Ayy = 0.0

	
	if swith_wall[1] == 1
		for j in 2:cellymax-1
			# 壁に並行な方向の速度を求めるためx-の壁面に対しvecAyを使用
			for i in 1:cellxmax
				tvAyx = 0.5*(vecAy[i,j-1,1]+vecAy[i,j,1])
				tvAyy = 0.5*(vecAy[i,j-1,2]+vecAy[i,j,2])
				Ayx = tvAyx / (tvAyx^2+tvAyy^2)^0.5
				Ayy = tvAyy / (tvAyx^2+tvAyy^2)^0.5
				nu = mu[i,j] / Qbase[i,j,1]
				ubar = abs(Ayx*Qbase[i,j,2] + Ayy*Qbase[i,j,3])
				yplus[i,j], up = cal_wall_val_spalding(ubar, wally[i,j], nu)
			end
		end
	end
	if swith_wall[2] == 1
		for j in 2:cellymax-1
			# 壁に並行な方向の速度を求めるためx+の壁面に対しvecAyを使用
			for i in 1:cellxmax
				tvAyx = 0.5*(vecAy[i,j-1,1]+vecAy[i,j,1])
				tvAyy = 0.5*(vecAy[i,j-1,2]+vecAy[i,j,2])
				Ayx = tvAyx / (tvAyx^2+tvAyy^2)^0.5
				Ayy = tvAyy / (tvAyx^2+tvAyy^2)^0.5
				nu = mu[i,j] / Qbase[i,j,1]
				ubar = abs(Ayx*Qbase[i,j,2] + Ayy*Qbase[i,j,3])
				yplus[i,j], up = cal_wall_val_spalding(ubar, wally[i,j], nu)
			end
		end
	end
	if swith_wall[3] == 1
		for i in 2:cellxmax-1
			for j in 1:cellymax
			# 壁に並行な方向の速度を求めるためy-の壁面に対しvecAxを使用
			
				tvAxx = 0.5*(vecAx[i-1,j,1]+vecAx[i,j,1])
				tvAxy = 0.5*(vecAx[i-1,j,2]+vecAx[i,j,2])
				Axx = tvAxx / (tvAxx^2+tvAxy^2)^0.5
				Axy = tvAxy / (tvAxx^2+tvAxy^2)^0.5
				nu = mu[i,j] / Qbase[i,j,1]
				ubar = abs(Axx*Qbase[i,j,2] + Axy*Qbase[i,j,3])
				yplus[i,j], up = cal_wall_val_spalding(ubar, wally[i,j], nu)

				if yplus[i,j] > 500
					break
				end
			end
		end
	end
	if swith_wall[4] == 1
		for j in 1:cellymax
			# 壁に並行な方向の速度を求めるためy-の壁面に対しvecAxを使用
			for i in 2:cellxmax-1
				tvAxx = 0.5*(vecAx[i-1,j,1]+vecAx[i,j,1])
				tvAxy = 0.5*(vecAx[i-1,j,2]+vecAx[i,j,2])
				Axx = tvAxx / (tvAxx^2+tvAxy^2)^0.5
				Axy = tvAxy / (tvAxx^2+tvAxy^2)^0.5
				nu = mu[i,j] / Qbase[i,j,1]
				ubar = abs(Axx*Qbase[i,j,2] + Axy*Qbase[i,j,3])
				yplus[i,j], up = cal_wall_val_spalding(ubar, wally[i,j], nu)
			end
		end
	end

	return yplus
end

function cal_wall_val_spalding(u,y,nu)
	# spalding 定数
	k = 0.4
	B = 5.5

	# ニュートン法
	up = 200.0
	for i in 1:100
		old = up
		
		up = up - spalding(up,u,y,nu,k,B)/spalding_dash(up,u,y,nu,k,B)

		if abs(up-old)/up < 10^(-4)
			break
		end
	end

	yplus = u*y/nu/up
	return yplus, up
end

function spalding(up,u,y,nu,k,B)
	t = k*up
	F = - 1/up* u*y/nu + up + exp(-k*B)*(exp(t)-1-t-t^2/2-t^3/6)
	return F
end

function spalding_dash(up,u,y,nu,k,B)
	t = k*up
	Fdash = 1/up^2* u*y/nu + 1 + exp(-k*B)*(k*exp(t)-k-k^2*up-k^3*up^2/2)
	return Fdash
end


# test
#=
uu = 6.64037328
yy = 0.06299976860
nunu = 1.406107133e-5
println(cal_wall_val_spalding(uu,yy,nunu))
throw(UndefVarError(:x))
=#

function wallf_Van_Driest(yplus)
	f = 1 - exp(-yplus/26)
	return f
end
