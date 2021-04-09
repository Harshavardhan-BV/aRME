using DataFrames, FreqTables, CSV, Base.Threads, CategoricalArrays

function onOffPair(base_rates::Array{Float64,1}, states::Array{Int,1}, lambda::Float64)
	snew = [0,0]
	for i in 1:2
		if i == 1
			s1 = states[1]
			s2 = states[2]
		else
			s1 = states[2]
			s2 = states[1]
		end

		if s1 == 1
			rates = base_rates[1:3]
			push!(rates,0)
		else
			rates = [0.0]
			append!(rates, base_rates[2:4])
		end
		if s2 == 1
			rates[2] *= lambda
			rates[1] *= lambda
		end
		rates = rates./sum(rates)
		rates = cumsum(rates)
		r = rand(1)[1]
		if r<rates[2]
			snew[i] = 1
		end
	end
	return snew
end

function prodDeg(RNAcount::Int, state::Int, pDeath::Float64)
	if state == 1
		RNAcount += 1
	end
	if RNAcount > 1 && rand(1)[1]<pDeath
		RNAcount -= 1
	end
	return RNAcount
end

function prodDegPair(RNAcounts::Array{Int, 1}, states::Array{Int, 1}, pDeath::Float64)
	RNAnew = [0,0]
	for i in 1:2
		RNAnew[i] = prodDeg(RNAcounts[i], states[i], pDeath)
	end
	return RNAnew
end

function simulationFullPair(base_rates::Array{Float64,1}, nSteps::Int, pDeath::Float64, 
	lambda::Float64)
	state_list = rand([0,1], 2)
	state_list = [state_list; zeros(Int, (nSteps-1),2)]
	RNAcount_list = zeros(Int, nSteps, 2)
	for i in 2:nSteps
		# RNAcount_list[i] = prodDegPair(RNAcount_list[i-1, :], state_list[i-1, :], pDeath)
		state_list[i] = onOffPair(base_rates, state_list[i-1, :], lambda)
	end
	return state_list, RNAcount_list
end

function simulationFinalPair(base_rates::Array{Float64,1}, nSteps::Int, pDeath::Float64,
	lambda::Float64)
	state = rand([0,1], 2)
	RNAcount = [0,0]
	for i in 2:nSteps
		RNAcount = prodDegPair(RNAcount, state, pDeath)
		state = onOffPair(base_rates, state,lambda)
	end
	return [state, RNAcount]
end



function simRepsPair(base_rates, nSteps::Int, nReps::Int, pDeath::Float64, full::Bool,
	lambda::Float64)
	if !full
		state_list = zeros(Int, nReps, 2)
		RNAcount_list = zeros(Int, nReps, 2)
		for i in 1:nReps
			x = simulationFinalPair(base_rates, nSteps, pDeath,lambda)
			state_list[i, :] = x[1]
			# RNAcount_list[i, :] = x[2]
		end
		df = DataFrame(p1 = levels!(CategoricalArray(state_list[:,1]), [0,1]),
			p2 = levels!(CategoricalArray(state_list[:,2]), [0,1]))
		p = freqtable(df, :p1, :p2).array./nReps
		p00 = p[1,1]
		p01 = p[1,2]
		p10 = p[2,1]
		p11 = p[2,2]
		return [nSteps, [p00, p01, p10, p11], state_list, RNAcount_list]
	else
		#### This part is off. Check this!!!
		#=state_list = zeros(Int, nSteps,2)
		RNAcount_list = zeros(Int, nSteps,2)
		for i in 1:nReps
			x = simulationFullPair(base_rates, nSteps, pDeath, lambda)
			state_list = [state_list x[1]]
			RNAcount_list = [RNAcount_list x[2]]
		end
		n = nReps+2
		state_list = state_list[:,3:n]
		RNAcount_list = RNAcount_list[:,3:n]
		plist = [freqtable(levels!(CategoricalArray(state_list[i,:]), [0,1])).array[2] for i in 1:nSteps]
		return [DataFrame(Steps = 1:nSteps, pOns = plist), state_list, 
		RNAcount_list]=#
		####
		println("Method under construction")
		return
	end
end

function baseRate()
	base_rates = rand(4)
	base_rates = base_rates./sum(base_rates)
	return base_rates
end

function multiRatesFinProbsDep(nRates::Int, nSteps::Int, nReps::Int,
	randDeath::Bool, csv::Bool=true, pDeath::Float64=0.5)
	df = DataFrame(stayOn = Float64[], on = Float64[], off = Float64[],
		stayOff = Float64[],  p00 = Float64[], p01 = Float64[],
		p10 = Float64[], p11 = Float64[], lambda=Float64[], pDeath = Float64[])

	baseRateList = [baseRate() for i in 1:nRates]
	if randDeath
		pDeathList = rand(nRates)
	else
		pDeathList = fill(pDeath, nRates)
	end
	lambdaList = 0.01 .+ 99.99*rand(nRates)
	Threads.@threads for i in 1:nRates
		x = simRepsPair(baseRateList[i], nSteps, nReps, pDeathList[i], false, 
			lambdaList[i])
		y = append!(baseRateList[i], x[2])
		y = append!(y, [lambdaList[i], pDeathList[i]])
		push!(df, y)
	end
	if csv
		CSV.write(join(["DependentProbs", nSteps, nReps, nRates, "csv"], "_", "."),
			df)
	end
	return df
end

chunk(arr, n) = [arr[i:min(i + n - 1, end)] for i in 1:n:length(arr)] 

function parallel(nRates::Int, nSteps::Int, nReps::Int,
	randDeath::Bool, ID::Int = 1,pDeath::Float64=0.5)
	n = nthreads()
	j = Int(round(nRates/n))
	chunks = chunk(1:nRates, j)
	dfList = []
	Threads.@threads for i in 1:n
		nR = length(chunks[i])
		println(nR)
		df = multiRatesFinProbsDep(nR, nSteps, nReps, randDeath, false, pDeath)
		push!(dfList, df)
	end
	df1 = dfList[1]
	for i in 2:n
		df1 = vcat(df1, dfList[i])
	end
	CSV.write(join(["DependentProbs", nSteps, nReps, nRates, ID, "csv"], "_", "."),
		df1)
	return df1
end

@time x = parallel(35000, 10000, 100, true)
