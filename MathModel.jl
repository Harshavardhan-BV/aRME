using DataFrames, FreqTables, CSV, Base.Threads, CategoricalArrays

function onOffPair(base_rates::Array{Float64,1}, states::Array{Int,1}, lambda::Float64)
	# Function to simulate Activation/Deactivation of gene pair
	snew = [0,0]
	for i in 1:2 
		if i == 1 # switch s1 and s2 depending on gene?
			s1 = states[1]
			s2 = states[2]
		else
			s1 = states[2]
			s2 = states[1]
		end

		if s1 == 1 # If said gene is active 
			rates = base_rates[1:3] 
			push!(rates,0) # If gene on make stayOff prob 0
		else
			rates = [0.0]
			append!(rates, base_rates[2:4]) # If gene off make stayOn prob 0
		end
		if s2 == 1 # If other gene is active scale on and stayOn rates by λ
			rates[2] *= lambda 
			rates[1] *= lambda 
		end
		rates = rates./sum(rates) # Normalize rates 
		rates = cumsum(rates) # Add on and stayOn probabilities
		r = rand(1)[1]
		if r<rates[2]
			snew[i] = 1 # Activate gene depending on RNG 
		end
	end
	return snew # Return new states
end

function prodDeg(RNAcount::Int, state::Int, pDeath::Float64)
	# Function to simulate RNA production/decay of a single gene
	if state == 1
		RNAcount += 1 # RNA increases if gene on
	end
	if RNAcount > 1 && rand(1)[1]<pDeath
		RNAcount -= 1 # RNA decreases by decay prob
	end
	return RNAcount # Return new RNA count
end

function prodDegPair(RNAcounts::Array{Int, 1}, states::Array{Int, 1}, pDeath::Float64)
	# Function to simulate RNA production/decay of gene pair
	RNAnew = [0,0] # Initialize RNA produced at zero count
	for i in 1:2
		RNAnew[i] = prodDeg(RNAcounts[i], states[i], pDeath) # Seperate function calls for each gene? 
	end
	return RNAnew # Return new RNA count
end

function simulationFullPair(base_rates::Array{Float64,1}, nSteps::Int, pDeath::Float64, 
	lambda::Float64)
	# Function to simulate gene pair with states saved at all timepoints (not active)
	state_list = rand([0,1], 2) # Initialise genes at random states
	state_list = [state_list; zeros(Int, (nSteps-1),2)] # Extend array for full time length
	RNAcount_list = zeros(Int, nSteps, 2) 
	for i in 2:nSteps
		# RNAcount_list[i] = prodDegPair(RNAcount_list[i-1, :], state_list[i-1, :], pDeath)
		state_list[i] = onOffPair(base_rates, state_list[i-1, :], lambda) # Call Gene switch simulation
	end
	return state_list, RNAcount_list # Return timeseries of state ~~and RNAcount~~
end

function simulationFinalPair(base_rates::Array{Float64,1}, nSteps::Int, pDeath::Float64,
	lambda::Float64)
	# Function to simulate gene pair with states saved at endpoint
	state = rand([0,1], 2) # Initialise genes at random states
	RNAcount = [0,0] # Initialize RNAs at zero count
	for i in 2:nSteps
		RNAcount = prodDegPair(RNAcount, state, pDeath) # Call RNA simulation
		state = onOffPair(base_rates, state,lambda) # Call Gene switch simulation
	end
	return [state, RNAcount] # Return only last state and RNAcount 
end

function simRepsPair(base_rates, nSteps::Int, nReps::Int, pDeath::Float64, full::Bool,
	lambda::Float64)
	# Function to simulate over replicates
	if !full
		state_list = zeros(Int, nReps, 2)
		RNAcount_list = zeros(Int, nReps, 2)
		for i in 1:nReps # Iterate over replicates
			x = simulationFinalPair(base_rates, nSteps, pDeath,lambda) # Call final pair RNA/Gene state simulation
			state_list[i, :] = x[1]
			# RNAcount_list[i, :] = x[2]
		end
		df = DataFrame(p1 = levels!(CategoricalArray(state_list[:,1]), [0,1]), # Categorical arrays with 0 and 1 as categories
			p2 = levels!(CategoricalArray(state_list[:,2]), [0,1])) # Split state_list by gene and add to dataframe
		p = freqtable(df, :p1, :p2).array./nReps # Construct frequency table from state list 
		p00 = p[1,1] # both inactive
		p01 = p[1,2] # 0 active 1 inactive
		p10 = p[2,1] # 0 inactive 1 active
		p11 = p[2,2] # both active
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
	# Function to generate rate parameters for on/off
	base_rates = rand(4)  #Random rates [stayOn,on,off,stayOff]
	base_rates = base_rates./sum(base_rates) # normalise rates
	return base_rates
end

function multiRatesFinProbsDep(nRates::Int, nSteps::Int, nReps::Int,
	randDeath::Bool, csv::Bool=true, pDeath::Float64=0.5)
	df = DataFrame(stayOn = Float64[], on = Float64[], off = Float64[],
		stayOff = Float64[],  p00 = Float64[], p01 = Float64[],
		p10 = Float64[], p11 = Float64[], lambda=Float64[], pDeath = Float64[]) # Initialise empty dataframe

	baseRateList = [baseRate() for i in 1:nRates] # Generate a list of switching rate
	if randDeath
		pDeathList = rand(nRates) # Generate a list of random death rates
	else
		pDeathList = fill(pDeath, nRates) # If not random fill everything with same probability
	end
	lambdaList = 0.01 .+ 99.99*rand(nRates) #lambda ??
	Threads.@threads for i in 1:nRates # Parallelize for loop
		x = simRepsPair(baseRateList[i], nSteps, nReps, pDeathList[i], false, 
			lambdaList[i])
		y = append!(baseRateList[i], x[2]) # Append rates and pij in list
		y = append!(y, [lambdaList[i], pDeathList[i]]) # Append lambda and probability of death
		push!(df, y) # Push the list to dataframe
	end
	if csv
		CSV.write(join(["DependentProbs", nSteps, nReps, nRates, "csv"], "_", "."),
			df) # Save dataframe as CSV if enabled
	end
	return df
end

chunk(arr, n) = [arr[i:min(i + n - 1, end)] for i in 1:n:length(arr)] # Function to split array to subarrays of size n (±1?)

function parallel(nRates::Int, nSteps::Int, nReps::Int,
	randDeath::Bool, ID::Int = 1,pDeath::Float64=0.5)
	n = nthreads() # Use all the threads
	j = Int(round(nRates/n)) 
	chunks = chunk(1:nRates, j) # Split rates by no of threads
	dfList = []
	Threads.@threads for i in 1:n # Parallelize for loop
		nR = length(chunks[i])
		println(nR)
		df = multiRatesFinProbsDep(nR, nSteps, nReps, randDeath, false, pDeath)
		push!(dfList, df) # Push the list to dataframe
	end
	df1 = dfList[1] 
	for i in 2:n
		df1 = vcat(df1, dfList[i]) # Parallelization concatanates horizontally - switching to vertical for column wise
	end
	CSV.write(join(["DependentProbs", nSteps, nReps, nRates, ID, "csv"], "_", "."), # Save as a CSV [name_format=TimeSteps-Replicates-NumberOfRates-ID]
		df1) # Save dataframe as CSV
	return df1
end

# Run simulation with nRates, nSteps, nReps, randDeath
@time x = parallel(35000, 10000, 100, true) 
