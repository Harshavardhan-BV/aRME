using DataFrames, Base.Threads, CSV

function switchGene(r::Matrix{Float64}, s::Vector{Bool}, λ::Vector{Float64})
    # Function to simulate Activation/Deactivation of gene pair
    λ = λ .* s[end:-1:1] # Use lambda if other gene is on else 0
    λ = hcat(λ,λ,zeros(2),zeros(2)) # Only change on parameters
    r = r .+ λ .* r # Add lambda factor to rates
    n1 , n2 = [3,4] - p .* 2 
    r = hcat(r[n1:4:end],r[n2:4:end])'# If gene is on choose stayOn & off else choose on & stayOff
    r = r ./ sum(r, dims=2) # Normalize over horizonatal sum
    RNG = rand(2)
    s .= (RNG .< r[:,1]) # Switch on gene if activation probability satisfied
    return s
end

function synRNA(RNA::Vector{Int}, s::Vector{Bool}, δ::Vector{Float64})
    # Function to simulate Production/Degradation of RNA
    RNA .+= s # Increment RNA by 1 if gene active
    RNG = rand(2)
    RNA .-= (RNG.<δ) .* (RNA.>1) # Decrease RNA by 1 if decay probability satisfied "and" RNA count higher than 1
    return RNA
end

function endPoint(tmax::Int, r::Matrix{Float64}, λ::Vector{Float64}, δ::Vector{Float64})
    #Function to simulate gene pair with states saved at endpoint
    s = rand(Vector{Bool}([0,1]), 2) # Initialise genes at random states
    RNA = zeros(Int64,2) # Initialize RNA at zero count
    for t in 2:tmax
        synRNA(RNA,s,δ) # Call RNA simulation
        switchGene(r, s, λ) # Call Gene switch simulation
    end
    return [s, RNA] # Return only last state and RNAcount 
end

function initParms(nRates, legacy::Bool=false, randDeath::Bool=true, pDeath::Float64=0.5)
    # Function to initialize rate parameters
    if legacy # legacy assumes rate parameters are symmetric over the pair
        d=1 # Only 1 set of values per gene pair 
    else 
        d=2 # Each set of values for each gene in pair [G1 G2]
    end
    r_list = rand(d,4,nRates) # Make a (dx4)xn random array for gene rates [G1, G2]? x [stayOn, on, off, stayOff] x n
    λ_list = 0.01 .+ 99.99 .* rand(d,nRates) # Make a (d)xn random array for influence parameter
    if randDeath
        δ_list = rand(d,nRates) # Make a (d)xn random array for RNA decay rates 
    else
        δ_list = fill(pDeath, d, nRates) # Make a (d)xn array with same probability for RNA decay rates
    end
    if legacy # For legacy duplicate values over the "vertical" dimension
        r_list = vcat(r_list,r_list) 
        λ_list = vcat(λ_list,λ_list) 
        δ_list = vcat(δ_list,δ_list)
    end
    return r_list, λ_list, δ_list
end

function overReps(nReps::Int, tmax::Int, r::Matrix{Float64}, λ::Vector{Float64}, δ::Vector{Float64})
    # Function to iterate over replicates
    rep_list = [] # Empty list for storing 
    for i in 1:nReps
        ep = endPoint(tmax, r, λ, δ)
        push!(rep_list, ep)
    end
    ## Convert to df???
    return df 
end

nRates = 35000
nSteps = 10000
nReps = 100
randDeath = true

