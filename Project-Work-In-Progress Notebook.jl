### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 5f7fb2f6-dbca-45ea-bb0a-1914260d0876
begin

	# load some external packages 
	using PlutoUI
	using DataFrames
	using BSON
	using GLPK
	using PrettyTables
	using Plots
	using CSV
	using Optim
	using Statistics
	using LinearAlgebra
	using Plots
	
	# setup my paths (where are my files?)
	_PATH_TO_ROOT = pwd() 
	_PATH_TO_SRC = joinpath(_PATH_TO_ROOT,"src")
	_PATH_TO_MODEL = joinpath(_PATH_TO_ROOT,"model")
	_PATH_TO_FIGS = joinpath(_PATH_TO_ROOT,"figs")
	
	# load the ENGRI 1120 project code library -
	include(joinpath(_PATH_TO_SRC,"Include.jl"))

	# load the model -
	MODEL = BSON.load(joinpath(_PATH_TO_MODEL,"model_v2.bson"), @__MODULE__)

	# show -
	nothing
end

# ╔═╡ 4855467c-3670-4e0b-a64e-0e09effa6e0d
md"""
## ENGRI 1120: Design and Analysis of a Sustainable Cell-Free Production Process for Industrially Important Small Molecules
"""

# ╔═╡ 8ca99e1f-afd6-4c1c-8918-db6d9747099c
html"""
<p style="font-size:20px;">Work In Progress: Matthew Chen, Elaina Varriale, Isabel Voellmicke, Avery Willis</br>
Smith School of Chemical and Biomolecular Engineering, Cornell University, Ithaca NY 14850</p>
"""

# ╔═╡ a0ad3474-1844-41bc-bd95-242aa94a5ff1
html"""
This design aims to optimize the production of 1,2 propanediol(PGDN) using the substrate maltose. PGDN has many uses. First, PGDN is used as a propellant in Otto fuel II, a type of fuel used to power torpedoes. PGDN is also used in medication to treat angina pectoris, a condition that includes chronic chest pain or discomfort. The substrate used to produce PGDN is maltose, a disaccharide sugar made of two glucose molecules. Maltose is purchased as dry crystals or syrup and is present in everyday foods such as cereals, certain fruits, and sweet potatoes. Potassium nitrate is also used as a reactant in this reaction. KNO3 is used in fertilizers, tree stump removal, rocket propellants and fireworks, and gunpowder. Finally, water is also used as a reactant. This design incorporates a 100µL reactor for $100 to house the reactions between reagents and products. These reactors have 2 input streams and 1 output stream. Additionally, “magical” separators for $20 will be used to split the desired products from byproducts produced in the reaction. This design produces a PGDN flow rate of 34.524 mmol/hr and has a total cost of $2880 and creates a revenue of $2,096,013.96 per year. 
"""

# ╔═╡ 40da982c-1cc4-4881-a2ea-fbeef5c46d2d
md"""
### Materials and Methods
"""

# ╔═╡ ae9106e9-2677-4596-9a97-1d0aa9f12152
md"""
##### Assumptions
"""

# ╔═╡ d643abdd-0467-4816-8b06-bf682e8f2404
html"""
<ul>
<li>Design Process will run 24/7</li>
<li>Cost to maintain system is negligible</li>
<li>No cost for shipping of materials</li>
<li>Process is isothermal and runs at an optimal temperature</li>
<li>Negligible power consumption</li>
<li>Water is free</li>
<li>Infinite volume syringes in syring pump</li>
<li>Ideal system at steady state</li>
</ul>
"""

# ╔═╡ 4867b51c-fb6c-42a9-b87d-4131e014b402
md"""
##### Setup the Flux Balance Analysis (FBA) calculation to estimate of optimal reaction rates
"""

# ╔═╡ 5f43cb39-e600-4530-9a3b-25e5251f78fd
html"""
Flux balance analysis is the mathematical approach to compute the needed flux for our reaction process. Our project team was given that the flux bounds were +/- 2.466 mmol per hour. With the given Flux Balance Analysis (FBA ) coding provided by Varner Labs, we input our given sugar feedstock and additional reactants and tested different inputs to find the optimal reaction rates. We determined that the maximum extent for our reaction was 2.466 mmol per hour. In addition, we were able to determine and check the input of Maltose, water, and KNO3 with the use of the stoichiometric coefficients of our reaction. Through the use of the stoichiometric coefficients and number of reactor chips, were we able to determine the input flow rates for our reactants to attain the needed amount of PDGN.
"""

# ╔═╡ 1a75041f-a7b5-491d-a74e-7e4740f0ceed
html"""
Futhermore, we set the flow rate for our syring pump to 0.01 mL/hr in order to approximate the amount of time that our input molecule spends in the reactor. Although we assume that our molecules will fully react, this was used to minimize the amount of water needed in the reaction. 
"""

# ╔═╡ 5432f738-c2cd-4727-821a-ca4fb4b04d19
begin

	# setup the FBA calculation for the project -

	# === SELECT YOUR PRODUCT HERE ==================================================== #
	# What rate are trying to maximize? (select your product)
	# rn:R08199 = isoprene
	# rn:28235c0c-ec00-4a11-8acb-510b0f2e2687 = PGDN
	# rn:rn:R09799 = Hydrazine
	# rn:R03119 = 3G
	idx_target_rate = find_reaction_index(MODEL,:reaction_number=>"rn:28235c0c-ec00-4a11-8acb-510b0f2e2687")
	# ================================================================================= #

	# First, let's build the stoichiometric matrix from the model object -
	(cia,ria,S) = build_stoichiometric_matrix(MODEL);

	# Next, what is the size of the system? (ℳ = number of metabolites, ℛ = number of reactions)
	(ℳ,ℛ) = size(S)

	# Next, setup a default bounds array => update specific elements
	# We'll correct the directionality below -
	Vₘ = (13.7)*(3600)*(50e-9)*(1000) # units: mmol/hr
	flux_bounds = [-Vₘ*ones(ℛ,1) Vₘ*ones(ℛ,1)]

	# update the flux bounds -> which fluxes can can backwards? 
	# do determine this: sgn(v) = -1*sgn(ΔG)
	updated_flux_bounds = update_flux_bounds_directionality(MODEL,flux_bounds)

	# hard code some bounds that we know -
	updated_flux_bounds[44,1] = 0.0  # ATP synthesis can't run backwards 

	# What is the default mol flow input array => update specific elements
	# strategy: start with nothing in both streams, add material(s) back
	n_dot_input_stream_1 = zeros(ℳ,1)	# stream 1
	n_dot_input_stream_2 = zeros(ℳ,1)	# stream 2

	# === YOU NEED TO CHANGE BELOW HERE ====================================================== #
	# Let's lookup stuff that we want/need to supply to the chip to get the reactiont to go -
	# what you feed *depends upon your product*
	compounds_that_we_need_to_supply_feed_1 = [
		"h2o", "maltose"
	]

	# what are the amounts that we need to supply to chip in feed stream 1 (units: mmol/hr)?
	mol_flow_values_feed_1 = [
		7Vₘ 	; # h2o mmol/hr
		7Vₘ 	; # maltose mmol/hr
	]

	# what is coming into feed stream 2?
	compounds_that_we_need_to_supply_feed_2 = [
		"potassium nitrate"
	]

	# let's always add Vₘ into feed stream 2
	mol_flow_values_feed_2 = [
		2Vₘ 		; # glycerol mmol/hr
	]
	
	
	# === YOU NEED TO CHANGE ABOVE HERE ====================================================== #

	# stream 1:
	idx_supply_stream_1 = Array{Int64,1}()
	for compound in compounds_that_we_need_to_supply_feed_1
		idx = find_compound_index(MODEL,:compound_name=>compound)
		push!(idx_supply_stream_1,idx)
	end

	# stream 2:
	idx_supply_stream_2 = Array{Int64,1}()
	for compound in compounds_that_we_need_to_supply_feed_2
		idx = find_compound_index(MODEL,:compound_name=>compound)
		push!(idx_supply_stream_2,idx)
	end
	
	# supply for stream 1 and stream 2
	n_dot_input_stream_1[idx_supply_stream_1] .= mol_flow_values_feed_1
	n_dot_input_stream_2[idx_supply_stream_2] .= mol_flow_values_feed_2
	
	# setup the species bounds array -
	species_bounds = [-1.0*(n_dot_input_stream_1.+n_dot_input_stream_2) 1000.0*ones(ℳ,1)]

	# Lastly, let's setup the objective function -
	c = zeros(ℛ)
	c[idx_target_rate] = -1.0

	# show -
	nothing
end

# ╔═╡ 2e275308-40d1-473a-9834-5df647b99e0a
md"""
###### Check: did calculation converge?
"""

# ╔═╡ c0d2722d-1b85-4bc0-841c-53a2a80a9aea
begin

	# compute the optimal flux -
	result = calculate_optimal_flux_distribution(S, updated_flux_bounds, species_bounds, c);

	# get the open extent vector -
	ϵ_dot = result.calculated_flux_array

	# what is the composition coming out of the first chip?
	n_dot_out_chip_1 = (n_dot_input_stream_1 + n_dot_input_stream_2 + S*ϵ_dot);

	# did this converge?
	with_terminal() do

		# get exit/status information from the solver -
		exit_flag = result.exit_flag
		status_flag = result.status_flag

		# display -
		println("Computed optimal flux distribution w/exit_flag = 0: $(exit_flag==0) and status_flag = 5: $(status_flag == 5)")
	end
end

# ╔═╡ c4b25914-554f-4870-a63b-86e05f6864bb
md"""
##### Step 2: Compute the output of chips i = 2, ..., N.
"""

# ╔═╡ 74c935ba-23b0-45a5-88c1-126b98b4cf06
begin

	# setup calculation for chips i = 2,....,N
	N = 14 # number of chips

	# initialize some space to store the mol flow rates -
	series_mol_state_array = zeros(ℳ,N)
	exit_flag_array = Array{Int64,1}()
	status_flag_array = Array{Int64,1}()

	# the initial col of this array is the output of from chip 1
	for species_index = 1:ℳ
		series_mol_state_array[species_index,1] = n_dot_out_chip_1[species_index]
	end
	
	# assumption: we *always* feed potassium nitrate into port 2 - so we only need to update the input flow into port 1
	for chip_index = 2:N

		# update the input into the chip -
		n_dot_input_port_1 = series_mol_state_array[:,chip_index - 1] 		# the input to chip j comes from j - 1
	
		# setup the species bounds array -
		species_bounds_next_chip = [-1.0*(n_dot_input_port_1.+n_dot_input_stream_2) 1000.0*ones(ℳ,1)]

		# run the optimal calculation -
		result_next_chip = calculate_optimal_flux_distribution(S, updated_flux_bounds, species_bounds_next_chip, c);

		# grab the status and exit flags ... so we can check all is right with the world ...
		push!(exit_flag_array, result_next_chip.exit_flag)
		push!(status_flag_array, result_next_chip.status_flag)

		# Get the flux from the result object -
		ϵ_dot_next_chip = result_next_chip.calculated_flux_array

		# compute the output from chip j = chip_index 
		n_dot_out_next_chip = (n_dot_input_port_1 + n_dot_input_stream_2 + S*ϵ_dot_next_chip);

		# copy this state vector into the state array 
		for species_index = 1:ℳ
			series_mol_state_array[species_index,chip_index] = n_dot_out_next_chip[species_index]
		end

		# go around again ...
	end
end

# ╔═╡ 07709469-b7a0-4c7b-9b92-1162efa14dd6
exit_flag_array

# ╔═╡ c37fa831-a359-425a-9a66-28bb589e1104
status_flag_array

# ╔═╡ e933ddd9-8fd8-416a-8710-a64d3eb36f79
md"""
###### Table 1: State table describing the exit composition (mmol/hr) for each chip. 

Each row of the table shows a different compound, while the columns show the mol flow rate for component i in the output from chip i. The last two columns show the mass flow rate and mass fraction for component i in the exit from chip N.
"""

# ╔═╡ 8a732899-7493-45da-bd6d-ecfba04f3ef1
begin

	# what chip r we looking at?
	n_dot_output_chip = series_mol_state_array[:,end]

	# get the array of MW -
	MW_array = MODEL[:compounds][!,:compound_mw]

	# convert the output mol stream to a mass stream -
	mass_dot_output = (n_dot_output_chip.*MW_array)*(1/1000)

	# what is the total coming out?
	total_mass_out = sum(mass_dot_output)
	
	# display code makes the table -
	with_terminal() do

		# what are the compound names and code strings? -> we can get these from the MODEL object 
		compound_name_strings = MODEL[:compounds][!,:compound_name]
		compound_id_strings = MODEL[:compounds][!,:compound_id]
		
		# how many molecules are in the state array?
		ℳ_local = length(compound_id_strings)
	
		# initialize some storage -
		number_of_cols = 3 + N + 2
		state_table = Array{Any,2}(undef,ℳ_local,number_of_cols)

		# get the uptake array from the result -
		uptake_array = result.uptake_array

		# populate the state table -
		for compound_index = 1:ℳ_local
			state_table[compound_index,1] = compound_index
			state_table[compound_index,2] = compound_name_strings[compound_index]
			state_table[compound_index,3] = compound_id_strings[compound_index]

			for chip_index = 1:N
				tmp_value = abs(series_mol_state_array[compound_index, chip_index])
				state_table[compound_index,chip_index + 3] = (tmp_value) <= 1e-6 ? 0.0 : 
					series_mol_state_array[compound_index, chip_index]
			end

			# show the mass -
			tmp_value = abs(mass_dot_output[compound_index])
			state_table[compound_index,(N + 3 + 1)] = (tmp_value) <= 1e-6 ? 0.0 : mass_dot_output[compound_index]

			# show the mass fraction -
			# show the mass -
			tmp_value = abs(mass_dot_output[compound_index])
			state_table[compound_index, (N + 3 + 2)] = (tmp_value) <= 1e-6 ? 0.0 : 	
				(1/total_mass_out)*mass_dot_output[compound_index]
		end

		# build the table header -
		id_header_row = Array{String,1}()
		units_header_row = Array{String,1}()

		# setup id row -
		push!(id_header_row, "i")
		push!(id_header_row, "name")
		push!(id_header_row, "id")
		for chip = 1:N
			push!(id_header_row, "Chip $(chip)")
		end
		push!(id_header_row, "m_dot")
		push!(id_header_row, "ωᵢ_output")

		# setup units header row -
		push!(units_header_row, "")
		push!(units_header_row, "")
		push!(units_header_row, "")
		for chip = 1:N
			push!(units_header_row, "mmol/hr")
		end
		push!(units_header_row, "g/hr")
		push!(units_header_row, "")
		
		# header row -
		state_table_header_row = (id_header_row, units_header_row)
		
		# write the table -
		pretty_table(state_table; header=state_table_header_row)
	end
end

# ╔═╡ 64daa21a-ac42-4b20-9e6b-ec2d19cd50fc
md"""
###### Table 2: Optimal reaction extent table computed by flux balance analysis. 

Each row corresponds to a reaction in the model where $\dot{\epsilon}_{i}$ denotes the optimal open reaction extent computed by flux balance analysis. 
"""

# ╔═╡ 4723ed10-8b65-46f6-b825-3e8bc43d6004
md"""
Reaction string format: $(@bind rxn_string_ver Select(["first"=>"KEGG", "second"=>"HUMAN"]))
"""

# ╔═╡ 8de4fb56-8d56-4251-ad5b-478bae38f727
with_terminal() do

	# initialize some storage -
	flux_table = Array{Any,2}(undef,ℛ,6)

	# what are the reaction strings? -> we can get these from the MODEL object 
	reaction_strings = MODEL[:reactions][!,:reaction_markup]
	reaction_id = MODEL[:reactions][!,:reaction_number]

	# translate the reaction string to HUMAN -
	human_rxn_string_array = translation_reaction_string_to_human(MODEL)

	# populate the state table -
	for reaction_index = 1:ℛ
		flux_table[reaction_index,1] = reaction_index
		flux_table[reaction_index,2] = reaction_id[reaction_index]
		
		if (rxn_string_ver == "first")
			flux_table[reaction_index,3] = reaction_strings[reaction_index]
		else
			flux_table[reaction_index,3] = human_rxn_string_array[reaction_index]
		end
		flux_table[reaction_index,4] = flux_bounds[reaction_index,1]
		flux_table[reaction_index,5] = flux_bounds[reaction_index,2]

		# clean up the display -
		tmp_value = abs(ϵ_dot[reaction_index])
		flux_table[reaction_index,6] = tmp_value < 1e-6 ? 0.0 : ϵ_dot[reaction_index]
	end

	# header row -
	flux_table_header_row = (["i","RID","R","ϵ₁_dot LB", "ϵ₁_dot UB", "ϵᵢ_dot"],
		["","","", "mmol/hr", "mmol/hr", "mmol/hr"]);
		
	# write the table -
	pretty_table(flux_table; header=flux_table_header_row)
end

# ╔═╡ 7720a5a6-5d7e-4c79-8aba-0a4bb04973af
md"""
##### Compute the downstream separation using Magical Separation Units (MSUs)
"""

# ╔═╡ 650b769e-6c70-43c4-8b3d-8e487bd66e3f
html"""
Given the desired output flow rate for our product, PDGN, which was about 1 g/hr with above 95% purity, we needed to find an optimal number of reactors and separators. Through testing various numbers within Julia notebook, we were able to attain the desired purities. This allowed us to sell our byproduces based on their purities"""

# ╔═╡ 78d9b6fe-53da-4944-bdb4-f88e767387c6
md"""
###### Product Separation Goals
PGDN - 99.7%
CO2 - 99.8%
KOH - 85%
Acetate - 95%
"""

# ╔═╡ 19977dec-e455-4069-a138-22b1ca8bf4d8
md"""
PGDN:
"""

# ╔═╡ 8fa0d539-5f6d-45c6-9151-f47273434f9c
# how many levels are we going to have in the separation tree?
number_of_levels = 7

# ╔═╡ 28c50656-7543-4ff2-b5cb-fd7810ec5e79
begin

	# define the split -
	θ = 0.75

	# most of the "stuff" has a 1 - θ in the up, and a θ in the down
	u = (1-θ)*ones(ℳ,1)
	d = θ*ones(ℳ,1)

	# However: the desired product has the opposite => correct for my compound of interest -> this is compound i = ⋆
	idx_target_compound = find_compound_index(MODEL,:compound_name=>"PGDN")

	# correct defaults -
	u[idx_target_compound] = θ
	d[idx_target_compound] = 1 - θ

	# let's compute the composition of the *always up* stream -
	
	# initialize some storage -
	species_mass_flow_array_top = zeros(ℳ,number_of_levels)
	species_mass_flow_array_bottom = zeros(ℳ,number_of_levels)

	for species_index = 1:ℳ
		value = mass_dot_output[species_index]
		species_mass_flow_array_top[species_index,1] = value
		species_mass_flow_array_bottom[species_index,1] = value
	end
	
	for level = 2:number_of_levels

		# compute the mass flows coming out of the top -
		m_dot_top = mass_dot_output.*(u.^(level-1))
		m_dot_bottom = mass_dot_output.*(d.^(level-1))

		# update my storage array -
		for species_index = 1:ℳ
			species_mass_flow_array_top[species_index,level] = m_dot_top[species_index]
			species_mass_flow_array_bottom[species_index,level] = m_dot_bottom[species_index]
		end
	end
	
	# what is the mass fraction in the top stream -
	species_mass_fraction_array_top = zeros(ℳ,number_of_levels)
	species_mass_fraction_array_bottom = zeros(ℳ,number_of_levels)

	# array to hold the *total* mass flow rate -
	total_mdot_top_array = zeros(number_of_levels)
	total_mdot_bottom_array = zeros(number_of_levels)
	
	# this is a dumb way to do this ... you're better than that JV come on ...
	T_top = sum(species_mass_flow_array_top,dims=1)
	T_bottom = sum(species_mass_flow_array_bottom,dims=1)
	for level = 1:number_of_levels

		# get the total for this level -
		T_level_top = T_top[level]
		T_level_bottom = T_bottom[level]

		# grab -
		total_mdot_top_array[level] = T_level_top
		total_mdot_bottom_array[level] = T_level_bottom

		for species_index = 1:ℳ
			species_mass_fraction_array_top[species_index,level] = (1/T_level_top)*
				(species_mass_flow_array_top[species_index,level])
			species_mass_fraction_array_bottom[species_index,level] = (1/T_level_bottom)*
				(species_mass_flow_array_bottom[species_index,level])
		end
	end
end

# ╔═╡ 7db73d0b-ea5f-4619-be61-a57676868be4
begin

	stages = (1:number_of_levels) |> collect
	plot(stages,species_mass_fraction_array_top[idx_target_compound,:], linetype=:steppre,lw=2,legend=:bottomright, 
		label="Mass fraction i = PDO Tops")
	xlabel!("Stage index l",fontsize=18)
	ylabel!("Tops mass fraction ωᵢ (dimensionless)",fontsize=18)

	# make a 0.997 line target line -
	target_line = 0.997*ones(number_of_levels)
	plot!(stages, target_line, color="red", lw=2,linestyle=:dash, label="Target 99.7% purity")
end

# ╔═╡ 519f62ab-75a9-44fb-bf01-560d0efa11e7
with_terminal() do

	# initialize some space -
	state_table = Array{Any,2}(undef, number_of_levels, 3)
	for level_index = 1:number_of_levels
		state_table[level_index,1] = level_index
		state_table[level_index,2] = species_mass_fraction_array_top[idx_target_compound, level_index]
		state_table[level_index,3] = total_mdot_top_array[level_index]
	end
	
	# header -
	state_table_header_row = (["stage","ωᵢ i=⋆ top","mdot"],
			["","","g/hr"]);

	# write the table -
	pretty_table(state_table; header=state_table_header_row)
end

# ╔═╡ b2088b5c-146d-4ad9-8384-c8ee05441362
md"""
CO2:
"""

# ╔═╡ d24b8ff2-6493-40e5-9431-5c7388b36ec5
begin
	number_of_levels2 = 9
	# most of the "stuff" has a 1 - θ in the up, and a θ in the down
	u2 = (1-θ)*ones(ℳ,1)
	d2 = θ*ones(ℳ,1)

	# However: the desired product has the opposite => correct for my compound of interest -> this is compound i = ⋆
	idx_target_compound2 = find_compound_index(MODEL,:compound_name=>"co2")

	# correct defaults -
	u2[idx_target_compound2] = θ
	d2[idx_target_compound2] = 1 - θ

	# let's compute the composition of the *always up* stream -
	
	# initialize some storage -
	species_mass_flow_array_top2 = zeros(ℳ,number_of_levels2)
	species_mass_flow_array_bottom2 = zeros(ℳ,number_of_levels2)
	mass_dot_output2 = mass_dot_output.*(d.^(1))

	for species_index = 1:ℳ
		value2 = mass_dot_output2[species_index]
		species_mass_flow_array_top2[species_index,1] = value2
		species_mass_flow_array_bottom2[species_index,1] = value2
	end
	
	for level = 2:number_of_levels2

		# compute the mass flows coming out of the top -
		m_dot_top2 = mass_dot_output2.*(u2.^(level-1))
		m_dot_bottom2 = mass_dot_output2.*(d2.^(level-1))

		# update my storage array -
		for species_index = 1:ℳ
			species_mass_flow_array_top2[species_index,level] = m_dot_top2[species_index]
			species_mass_flow_array_bottom2[species_index,level] = m_dot_bottom2[species_index]
		end
	end
	
	# what is the mass fraction in the top stream -
	species_mass_fraction_array_top2 = zeros(ℳ,number_of_levels2)
	species_mass_fraction_array_bottom2 = zeros(ℳ,number_of_levels2)

	# array to hold the *total* mass flow rate -
	total_mdot_top_array2 = zeros(number_of_levels2)
	total_mdot_bottom_array2 = zeros(number_of_levels2)
	
	# this is a dumb way to do this ... you're better than that JV come on ...
	T_top2 = sum(species_mass_flow_array_top2,dims=1)
	T_bottom2 = sum(species_mass_flow_array_bottom2,dims=1)
	for level = 1:number_of_levels2

		# get the total for this level -
		T_level_top2 = T_top2[level]
		T_level_bottom2 = T_bottom2[level]

		# grab -
		total_mdot_top_array2[level] = T_level_top2
		total_mdot_bottom_array2[level] = T_level_bottom2

		for species_index = 1:ℳ
			species_mass_fraction_array_top2[species_index,level] = (1/T_level_top2)*
				(species_mass_flow_array_top2[species_index,level])
			species_mass_fraction_array_bottom2[species_index,level] = (1/T_level_bottom2)*
				(species_mass_flow_array_bottom2[species_index,level])
		end
	end
end

# ╔═╡ 379c70d8-f4a4-4603-9aab-b3ab20f92aaa
begin

	stages2 = (1:number_of_levels2) |> collect
	plot(stages2,species_mass_fraction_array_top2[idx_target_compound2,:], linetype=:steppre,lw=2,legend=:bottomright, 
		label="Mass fraction i = PDO Tops")
	xlabel!("Stage index l - 1",fontsize=18)
	ylabel!("Tops mass fraction ωᵢ (dimensionless)",fontsize=18)

	# make a 0.998 line target line -
	target_line2 = 0.998*ones(number_of_levels2)
	plot!(stages2, target_line2, color="red", lw=2,linestyle=:dash, label="Target 99.8% purity")
end

# ╔═╡ d8ddf8ac-82ef-49f2-b512-0ec4fe79aa4b
md"""
Table starts at stage 2 because stage 1 was the downstream separator flow from the initial PGDN separator
"""

# ╔═╡ 5dcd4195-259c-4008-adbb-b4c4abebba1c
with_terminal() do

	# initialize some space -
	state_table2 = Array{Any,2}(undef, number_of_levels2, 3)
	for level_index = 1:number_of_levels2
		state_table2[level_index,1] = level_index + 1
		state_table2[level_index,2] = species_mass_fraction_array_top2[idx_target_compound2, level_index]
		state_table2[level_index,3] = total_mdot_top_array2[level_index]
	end
	
	# header -
	state_table_header_row = (["stage","ωᵢ i=⋆ top","mdot"],
			["","","g/hr"]);

	# write the table -
	pretty_table(state_table2; header=state_table_header_row)
end

# ╔═╡ 6bb624af-d117-416d-b38c-750cd192c182
md"""
KOH:
"""

# ╔═╡ a675e809-c596-4783-9f4a-0a6e1cf41cd2
begin
	number_of_levels3 = 3
	# most of the "stuff" has a 1 - θ in the up, and a θ in the down
	u3 = (1-θ)*ones(ℳ,1)
	d3 = θ*ones(ℳ,1)

	# However: the desired product has the opposite => correct for my compound of interest -> this is compound i = ⋆
	idx_target_compound3 = find_compound_index(MODEL,:compound_name=>"potassium hydroxide")

	# correct defaults -
	u3[idx_target_compound3] = θ
	d3[idx_target_compound3] = 1 - θ

	# let's compute the composition of the *always up* stream -
	
	# initialize some storage -
	species_mass_flow_array_top3 = zeros(ℳ,number_of_levels3)
	species_mass_flow_array_bottom3 = zeros(ℳ,number_of_levels3)
	mass_dot_output3 = mass_dot_output2.*(d2.^(1))

	for species_index = 1:ℳ
		value3 = mass_dot_output3[species_index]
		species_mass_flow_array_top3[species_index,1] = value3
		species_mass_flow_array_bottom3[species_index,1] = value3
	end
	
	for level = 2:number_of_levels3

		# compute the mass flows coming out of the top -
		m_dot_top3 = mass_dot_output3.*(u3.^(level-1))
		m_dot_bottom3 = mass_dot_output3.*(d3.^(level-1))

		# update my storage array -
		for species_index = 1:ℳ
			species_mass_flow_array_top3[species_index,level] = m_dot_top3[species_index]
			species_mass_flow_array_bottom3[species_index,level] = m_dot_bottom3[species_index]
		end
	end
	
	# what is the mass fraction in the top stream -
	species_mass_fraction_array_top3 = zeros(ℳ,number_of_levels3)
	species_mass_fraction_array_bottom3 = zeros(ℳ,number_of_levels3)

	# array to hold the *total* mass flow rate -
	total_mdot_top_array3 = zeros(number_of_levels3)
	total_mdot_bottom_array3 = zeros(number_of_levels3)
	
	# this is a dumb way to do this ... you're better than that JV come on ...
	T_top3 = sum(species_mass_flow_array_top3,dims=1)
	T_bottom3 = sum(species_mass_flow_array_bottom3,dims=1)
	for level = 1:number_of_levels3

		# get the total for this level -
		T_level_top3 = T_top3[level]
		T_level_bottom3 = T_bottom3[level]

		# grab -
		total_mdot_top_array3[level] = T_level_top3
		total_mdot_bottom_array3[level] = T_level_bottom3

		for species_index = 1:ℳ
			species_mass_fraction_array_top3[species_index,level] = (1/T_level_top3)*
				(species_mass_flow_array_top3[species_index,level])
			species_mass_fraction_array_bottom3[species_index,level] = (1/T_level_bottom3)*
				(species_mass_flow_array_bottom3[species_index,level])
		end
	end
end

# ╔═╡ 4f92b0d0-bde0-47f4-80a6-376c3225b44b
begin

	stages3 = (1:number_of_levels3) |> collect
	plot(stages3,species_mass_fraction_array_top3[idx_target_compound3,:], linetype=:steppre,lw=2,legend=:bottomright, 
		label="Mass fraction i = PDO Tops")
	xlabel!("Stage index l - 2",fontsize=18)
	ylabel!("Tops mass fraction ωᵢ (dimensionless)",fontsize=18)

	# make a 0.85 line target line -
	target_line3 = 0.85*ones(number_of_levels3)
	plot!(stages3, target_line3, color="red", lw=2,linestyle=:dash, label="Target 85% purity")
end

# ╔═╡ 56226d6f-56a3-45c7-9940-99c86c1140eb
with_terminal() do

	# initialize some space -
	state_table3 = Array{Any,2}(undef, number_of_levels3, 3)
	for level_index = 1:number_of_levels3
		state_table3[level_index,1] = level_index + 2
		state_table3[level_index,2] = species_mass_fraction_array_top3[idx_target_compound3, level_index]
		state_table3[level_index,3] = total_mdot_top_array3[level_index]
	end
	
	# header -
	state_table_header_row = (["stage","ωᵢ i=⋆ top","mdot"],
			["","","g/hr"]);

	# write the table -
	pretty_table(state_table3; header=state_table_header_row)
end

# ╔═╡ 0d29f2f0-3045-4810-83d2-35464a7ce6f3
begin
	number_of_levels4 = 5
	# most of the "stuff" has a 1 - θ in the up, and a θ in the down
	u4 = (1-θ)*ones(ℳ,1)
	d4 = θ*ones(ℳ,1)

	# However: the desired product has the opposite => correct for my compound of interest -> this is compound i = ⋆
	idx_target_compound4 = find_compound_index(MODEL,:compound_name=>"acetate")

	# correct defaults -
	u4[idx_target_compound4] = θ
	d4[idx_target_compound4] = 1 - θ

	# let's compute the composition of the *always up* stream -
	
	# initialize some storage -
	species_mass_flow_array_top4 = zeros(ℳ,number_of_levels4)
	species_mass_flow_array_bottom4 = zeros(ℳ,number_of_levels4)
	mass_dot_output4 = mass_dot_output3.*(d3.^(1))

	for species_index = 1:ℳ
		value4 = mass_dot_output4[species_index]
		species_mass_flow_array_top4[species_index,1] = value4
		species_mass_flow_array_bottom4[species_index,1] = value4
	end
	
	for level = 2:number_of_levels4

		# compute the mass flows coming out of the top -
		m_dot_top4 = mass_dot_output4.*(u4.^(level-1))
		m_dot_bottom4 = mass_dot_output4.*(d4.^(level-1))

		# update my storage array -
		for species_index = 1:ℳ
			species_mass_flow_array_top4[species_index,level] = m_dot_top4[species_index]
			species_mass_flow_array_bottom4[species_index,level] = m_dot_bottom4[species_index]
		end
	end
	
	# what is the mass fraction in the top stream -
	species_mass_fraction_array_top4 = zeros(ℳ,number_of_levels4)
	species_mass_fraction_array_bottom4 = zeros(ℳ,number_of_levels4)

	# array to hold the *total* mass flow rate -
	total_mdot_top_array4 = zeros(number_of_levels4)
	total_mdot_bottom_array4 = zeros(number_of_levels4)
	
	# this is a dumb way to do this ... you're better than that JV come on ...
	T_top4 = sum(species_mass_flow_array_top4,dims=1)
	T_bottom4 = sum(species_mass_flow_array_bottom4,dims=1)
	for level = 1:number_of_levels4

		# get the total for this level -
		T_level_top4 = T_top4[level]
		T_level_bottom4 = T_bottom4[level]

		# grab -
		total_mdot_top_array4[level] = T_level_top4
		total_mdot_bottom_array4[level] = T_level_bottom4

		for species_index = 1:ℳ
			species_mass_fraction_array_top4[species_index,level] = (1/T_level_top4)*
				(species_mass_flow_array_top4[species_index,level])
			species_mass_fraction_array_bottom4[species_index,level] = (1/T_level_bottom4)*
				(species_mass_flow_array_bottom4[species_index,level])
		end
	end
end

# ╔═╡ b8749f56-7038-4c7c-b6b5-ea4d05486ac2
begin

	stages4 = (1:number_of_levels4) |> collect
	plot(stages4,species_mass_fraction_array_top4[idx_target_compound4,:], linetype=:steppre,lw=2,legend=:bottomright, 
		label="Mass fraction i = PDO Tops")
	xlabel!("Stage index l - 3",fontsize=18)
	ylabel!("Tops mass fraction ωᵢ (dimensionless)",fontsize=18)

	# make a 0.95 line target line -
	target_line4 = 0.95*ones(number_of_levels4)
	plot!(stages4, target_line4, color="red", lw=2,linestyle=:dash, label="Target 95% purity")
end

# ╔═╡ cf0dc9cf-77ed-4304-8e0d-be6d2fd69fdd
with_terminal() do

	# initialize some space -
	state_table4 = Array{Any,2}(undef, number_of_levels4, 3)
	for level_index = 1:number_of_levels4
		state_table4[level_index,1] = level_index + 3
		state_table4[level_index,2] = species_mass_fraction_array_top4[idx_target_compound4, level_index]
		state_table4[level_index,3] = total_mdot_top_array4[level_index]
	end
	
	# header -
	state_table_header_row = (["stage","ωᵢ i=⋆ top","mdot"],
			["","","g/hr"]);

	# write the table -
	pretty_table(state_table4; header=state_table_header_row)
end

# ╔═╡ f58c8483-1c6f-44db-92c1-26dc8694e0fe
md"""
#### Net Present Value (NPV):
Compared to a fixed rate 10% per year investment
"""

# ╔═╡ 8fad738d-18fa-43c8-87ee-773859d6b07b
begin
	T = 10
	discount_rate = 0.1
	CF_array = [-2.880, 2096.01396, 2096.01396, 2096.01396, 2096.01396, 2096.01396, 2096.01396, 2096.01396, 2096.01396, 2096.01396, 2096.01396]
end

# ╔═╡ 954c7d60-a17d-4ce3-b4e7-0ae71a7b8b84
begin
	
	# main -
	discount_array = ones(T+1)
	for time_index = 2:(T+1)
    	tmp_term = (1+discount_rate)^(time_index-1)
    	discount_array[time_index] = tmp_term
	end
	
	# return -
	nothing
end

# ╔═╡ afceae22-4433-4fa1-8cc3-f876ef9a338a
# what is the PV?
PV = CF_array.*(1.0./discount_array)

# ╔═╡ 4506e0f0-b019-411b-9701-99bea11edfd4
# what is the NPV?
NPV = sum(PV)

# ╔═╡ 7a2ddad5-0162-45af-ba89-1028e0405144
md"""
### Results and Dicussion
"""

# ╔═╡ 77107b01-dae9-4900-b9f7-f0c0224a492b
html"""
Our team had the decision whether to run our microfluidic chips in parallel--in which equal feeds are fed into identical chips that operate in parallel and then recombined-- or in series-- in which the output of one reactor is the input of another. Using the provided Julia notebooks and coding for chips in parallel and in series, we determined that there were no significant changes in the output flow rates between both options.
"""

# ╔═╡ 7252cc6f-d4a7-42e6-a998-fdb0e91aa9b0
md"""
##### Chips in parallel:
"""

# ╔═╡ ec3afa97-ecaf-4370-b0c8-2b8e321e9e25
begin

	# setup the FBA calculation for the project -

	# === SELECT YOUR PRODUCT HERE ==================================================== #
	# What rate are trying to maximize? (select your product)
	# rn:R08199 = isoprene
	# rn:28235c0c-ec00-4a11-8acb-510b0f2e2687 = PGDN
	# rn:rn:R09799 = Hydrazine
	# rn:R03119 = 3G
	idx_target_rate2 = find_reaction_index(MODEL,:reaction_number=>"rn:28235c0c-ec00-4a11-8acb-510b0f2e2687")
	# ================================================================================= #

	# First, let's build the stoichiometric matrix from the model object -
	(cia2,ria2,S2) = build_stoichiometric_matrix(MODEL);

	# Next, what is the size of the system? (ℳ = number of metabolites, ℛ = number of reactions)
	(ℳ2,ℛ2) = size(S2)

	# Next, setup a default bounds array => update specific elements
	# We'll correct the directionality below -
	Vₘ2 = (13.7)*(3600)*(50e-9)*(1000) # units: mmol/hr
	flux_bounds2 = [-Vₘ2*ones(ℛ2,1) Vₘ2*ones(ℛ2,1)]

	# update the flux bounds -> which fluxes can can backwards? 
	# do determine this: sgn(v) = -1*sgn(ΔG)
	updated_flux_bounds2 = update_flux_bounds_directionality(MODEL,flux_bounds2)

	# hard code some bounds that we know -
	updated_flux_bounds2[44,1] = 0.0  # ATP synthesis can't run backwards 

	# What is the default mol flow input array => update specific elements
	# strategy: start with nothing in both streams, add material(s) back
	n_dot_input_stream_12 = zeros(ℳ2,1)	# stream 1
	n_dot_input_stream_22 = zeros(ℳ2,1)	# stream 2

	# === YOU NEED TO CHANGE BELOW HERE ====================================================== #
	# Let's lookup stuff that we want/need to supply to the chip to get the reactiont to go -
	# what you feed *depends upon your product*
	compounds_that_we_need_to_supply_feed_12 = [
		"h2o", "maltose"
	]

	# what are the amounts that we need to supply to chip in feed stream 1 (units: mmol/hr)?
	mol_flow_values_feed_12 = [
		Vₘ2/2 	; # h2o mmol/hr
		Vₘ2/2 	; # maltose mmol/hr (maybe: 0.822 or 6.1?)
	]

	# what is coming into feed stream 2?
	compounds_that_we_need_to_supply_feed_22 = [
		"potassium nitrate"
	]

	# let's always add Vₘ into feed stream 2
	mol_flow_values_feed_22 = [
		2Vₘ2 		; # potassium nitrate mmol/hr
	]
	
	
	# === YOU NEED TO CHANGE ABOVE HERE ====================================================== #

	# stream 1:
	idx_supply_stream_12 = Array{Int64,1}()
	for compound in compounds_that_we_need_to_supply_feed_12
		idx = find_compound_index(MODEL,:compound_name=>compound)
		push!(idx_supply_stream_12,idx)
	end

	# stream 2:
	idx_supply_stream_22 = Array{Int64,1}()
	for compound in compounds_that_we_need_to_supply_feed_22
		idx = find_compound_index(MODEL,:compound_name=>compound)
		push!(idx_supply_stream_22,idx)
	end
	
	# supply for stream 1 and stream 2
	n_dot_input_stream_12[idx_supply_stream_12] .= mol_flow_values_feed_12
	n_dot_input_stream_22[idx_supply_stream_22] .= mol_flow_values_feed_22
	
	# setup the species bounds array -
	species_bounds2 = [-1.0*(n_dot_input_stream_12.+n_dot_input_stream_22) 1000.0*ones(ℳ2,1)]

	# Lastly, let's setup the objective function -
	c2 = zeros(ℛ2)
	c2[idx_target_rate2] = -1.0

	# show -
	nothing
end

# ╔═╡ 19b28401-ee7a-4d9a-bafd-493c4cf70e84
begin

	# compute the optimal flux -
	result2 = calculate_optimal_flux_distribution(S2, updated_flux_bounds2, species_bounds2, c2);

	# get the open extent vector -
	ϵ_dot2 = result2.calculated_flux_array

	# what is the composition coming out of the first chip?
	n_dot_out_chip_12 = (n_dot_input_stream_12 + n_dot_input_stream_22 + S2*ϵ_dot2);

	# did this converge?
	with_terminal() do

		# get exit/status information from the solver -
		exit_flag = result2.exit_flag
		status_flag = result2.status_flag

		# display -
		println("Computed optimal flux distribution w/exit_flag = 0: $(exit_flag==0) and status_flag = 5: $(status_flag == 5)")
	end
end

# ╔═╡ 95c0b4f2-a01c-47b7-bb6e-64fbea5de550
md"""
###### Table 1: State table from a single chip (species mol flow rate mmol/hr at exit)
"""

# ╔═╡ 42fb46e3-3c37-483d-bafc-f4b9f4b11d1d
md"""
##### Financials
"""

# ╔═╡ b63d7fe7-845f-4f70-a5f8-81cfa5e1c018
html"""
Cost Per Hour for Materials
<ul>
<li>Maltose: $0.74</li>
<li>KNO3: $0.44</li>
<li>PGDN: $79.10</li>
<li>KOH: $1.45</li>
<li>Acetate: $159.74</li>
<li>Carbon Dioxide: $16.10</li>
</ul>

Fixed Cost: $2880
<ul>
<li>14 Chips at $100 each</li>
<li>Syring pump at $1000</li>
<li>24 Separators at $20</li>
</ul>

Hourly Cost: $1.18
Hourly Revenue: $240.45 
Hourly Profit: $239.27
"""

# ╔═╡ 31eefb26-ac25-41f3-92ef-9a7600d6a093
html"""
Annual Profit: $2,096,013.96
"""

# ╔═╡ ab9084e1-57a2-4c6b-ad32-8fee8c142c43
md"""
### Conclusions
"""

# ╔═╡ e6e1b8fd-05f5-4f86-a1d1-8a93f7a7854c
html"""
Overall, this design is favorable: it produces revenue that much outweighs the cost invested in the product. Additionally, the yield is slightly above the goal target of 1.0g/hr and meets the 95% purity standard. When examining our Net Present Value analysis for our initial investment of $2880, it is evident that over 10 years our design process is financially favorable when compared to a 10% annual investment rate.
"""

# ╔═╡ 37b61b8c-d983-4d67-938c-cd17b9749426
html"""
Relevant Links 
<ul>
<li>“Acetic acid 96%.” Sigma Aldrich, https://www.sigmaaldrich.com/US/en/product/mm/100062?context=product. Accessed 9 December 2021.</li>
<li>Agency for Toxic Substances and Disease Registry. Otto Fuel II and Its Components. September 1996, https://www.atsdr.cdc.gov/toxfaqs/tfacts77.pdf.</li>
<li>“Angina Pectoris.” Johns Hopkins Medicine, https://www.hopkinsmedicine.org/health/conditions-and-diseases/angina-pectoris. Accessed 9 December 2021.</li>
<li>“Carbon Dioxide.” Sigma-Aldrich, https://www.sigmaaldrich.com/US/en/search/co2?focus=products&page=1&perPage=30&sort=relevance&term=co2&type=product. Accessed 9 December 2021.</li>
<li>“Maltose, Reagent Grade, 100 g.” Carolina, https://www.carolina.com/catalog/detail.jsp?prodId=873750&utm_source=google&utm_medium=cpc&scid=scplp873750&sc_intid=873750&gclid=Cj0KCQiAqbyNBhC2ARIsALDwAsAeSxtll8DNKDJvkOxmh7yUEJdVwpXh_laVlKm9pGXHSSIsS3V-BH4aAhlOEALw_wcB. Accessed 9 December 2021.</li>
<li>National Center for Biotechnology Information. “PubChem Compound Summary for CID 22933, Propylene glycol dinitrate.” PubChem, https://pubchem.ncbi.nlm.nih.gov/compound/Propylene-glycol-dinitrate#section=Mechanism-of-Action. Accessed 9 December 2021.</li>
<li>“1,2-Propanediol analytical standard 57-55-6.” Sigma-Aldrich, https://www.sigmaaldrich.com/US/en/product/sial/12279. Accessed 9 December 2021.</li>
<li>“Potassium Hydroxide Solution.” Sigma-Aldrich, https://www.sigmaaldrich.com/US/en/search/koh?focus=products&page=1&perPage=30&sort=relevance&term=koh&type=product. Accessed 9 December 2021.</li>
<li>“Potassium nitrate BioReagent, suitable for cell culture, suitable for plant cell culture | 7757-79-1.” Sigma-Aldrich, https://www.sigmaaldrich.com/US/en/product/sigma/p8291. Accessed 9 December 2021.</li> 
<li>“Potassium Nitrate Pharmaceutical Secondary Standard; Certified Reference Material | 7757-79-1.” Sigma-Aldrich, https://www.sigmaaldrich.com/US/en/product/sial/phr2203. Accessed 9 December 2021.</li>
</ul>

Prices: 
<li>Maltose - https://www.carolina.com/catalog/detail.jsp?prodId=873750&utm_source=google&utm_medium=cpc&scid=scplp873750&sc_intid=873750&gclid=Cj0KCQiAqbyNBhC2ARIsALDwAsAeSxtll8DNKDJvkOxmh7yUEJdVwpXh_laVlKm9pGXHSSIsS3V-BH4aAhlOEALw_wcB</li>
<li>Potassium Nitrate - https://www.sigmaaldrich.com/US/en/product/sigma/p8291, used price for 25kg of this</li>
<li>Carbon Dioxide - https://www.sigmaaldrich.com/US/en/product/aldrich/295108<</li>
<li>Potassium Hydroxide - https://www.sigmaaldrich.com/US/en/product/sigald/221473, used price for 25g of this</li>
<li>Acetate - https://vitasmlab.biz/store, used middle price for 500mg of acetate on this website
</li>
<li>PGDN - approximated using 1mL of https://www.sigmaaldrich.com/US/en/product/sial/12279 and https://www.sigmaaldrich.com/US/en/product/sial/phr2203</li>
"""

# ╔═╡ 18b29a1a-4787-11ec-25e3-5f29ebd21430
html"""
<style>

main {
    max-width: 1200px;
    width: 85%;
    margin: auto;
    font-family: "Roboto, monospace";
}

a {
    color: blue;
    text-decoration: none;
}

.H1 {
    padding: 0px 30px;
}
</style"""

# ╔═╡ 16dca67c-f280-4a6f-bd79-308cf63dabf6
html"""
<script>

	// initialize -
	var section = 0;
	var subsection = 0;
	var headers = document.querySelectorAll('h3, h4');
	
	// main loop -
	for (var i=0; i < headers.length; i++) {
	    
		var header = headers[i];
	    var text = header.innerText;
	    var original = header.getAttribute("text-original");
	    if (original === null) {
	        
			// Save original header text
	        header.setAttribute("text-original", text);
	    } else {
	        
			// Replace with original text before adding section number
	        text = header.getAttribute("text-original");
	    }
	
	    var numbering = "";
	    switch (header.tagName) {
	        case 'H3':
	            section += 1;
	            numbering = section + ".";
	            subsection = 0;
	            break;
	        case 'H4':
	            subsection += 1;
	            numbering = section + "." + subsection;
	            break;
	    }

		// update the header text 
		header.innerText = numbering + " " + text;
	};
</script>
"""

# ╔═╡ Cell order:
# ╟─4855467c-3670-4e0b-a64e-0e09effa6e0d
# ╟─8ca99e1f-afd6-4c1c-8918-db6d9747099c
# ╟─a0ad3474-1844-41bc-bd95-242aa94a5ff1
# ╟─40da982c-1cc4-4881-a2ea-fbeef5c46d2d
# ╟─ae9106e9-2677-4596-9a97-1d0aa9f12152
# ╟─d643abdd-0467-4816-8b06-bf682e8f2404
# ╟─4867b51c-fb6c-42a9-b87d-4131e014b402
# ╟─5f43cb39-e600-4530-9a3b-25e5251f78fd
# ╟─1a75041f-a7b5-491d-a74e-7e4740f0ceed
# ╠═5432f738-c2cd-4727-821a-ca4fb4b04d19
# ╟─2e275308-40d1-473a-9834-5df647b99e0a
# ╠═c0d2722d-1b85-4bc0-841c-53a2a80a9aea
# ╠═c4b25914-554f-4870-a63b-86e05f6864bb
# ╠═74c935ba-23b0-45a5-88c1-126b98b4cf06
# ╠═07709469-b7a0-4c7b-9b92-1162efa14dd6
# ╠═c37fa831-a359-425a-9a66-28bb589e1104
# ╟─e933ddd9-8fd8-416a-8710-a64d3eb36f79
# ╠═8a732899-7493-45da-bd6d-ecfba04f3ef1
# ╟─64daa21a-ac42-4b20-9e6b-ec2d19cd50fc
# ╟─4723ed10-8b65-46f6-b825-3e8bc43d6004
# ╟─8de4fb56-8d56-4251-ad5b-478bae38f727
# ╟─7720a5a6-5d7e-4c79-8aba-0a4bb04973af
# ╟─650b769e-6c70-43c4-8b3d-8e487bd66e3f
# ╠═78d9b6fe-53da-4944-bdb4-f88e767387c6
# ╟─19977dec-e455-4069-a138-22b1ca8bf4d8
# ╟─8fa0d539-5f6d-45c6-9151-f47273434f9c
# ╠═28c50656-7543-4ff2-b5cb-fd7810ec5e79
# ╠═7db73d0b-ea5f-4619-be61-a57676868be4
# ╠═519f62ab-75a9-44fb-bf01-560d0efa11e7
# ╠═b2088b5c-146d-4ad9-8384-c8ee05441362
# ╠═d24b8ff2-6493-40e5-9431-5c7388b36ec5
# ╠═379c70d8-f4a4-4603-9aab-b3ab20f92aaa
# ╟─d8ddf8ac-82ef-49f2-b512-0ec4fe79aa4b
# ╠═5dcd4195-259c-4008-adbb-b4c4abebba1c
# ╟─6bb624af-d117-416d-b38c-750cd192c182
# ╠═a675e809-c596-4783-9f4a-0a6e1cf41cd2
# ╠═4f92b0d0-bde0-47f4-80a6-376c3225b44b
# ╠═56226d6f-56a3-45c7-9940-99c86c1140eb
# ╠═0d29f2f0-3045-4810-83d2-35464a7ce6f3
# ╠═b8749f56-7038-4c7c-b6b5-ea4d05486ac2
# ╠═cf0dc9cf-77ed-4304-8e0d-be6d2fd69fdd
# ╟─f58c8483-1c6f-44db-92c1-26dc8694e0fe
# ╠═8fad738d-18fa-43c8-87ee-773859d6b07b
# ╠═954c7d60-a17d-4ce3-b4e7-0ae71a7b8b84
# ╠═afceae22-4433-4fa1-8cc3-f876ef9a338a
# ╠═4506e0f0-b019-411b-9701-99bea11edfd4
# ╠═7a2ddad5-0162-45af-ba89-1028e0405144
# ╟─77107b01-dae9-4900-b9f7-f0c0224a492b
# ╠═7252cc6f-d4a7-42e6-a998-fdb0e91aa9b0
# ╠═ec3afa97-ecaf-4370-b0c8-2b8e321e9e25
# ╠═19b28401-ee7a-4d9a-bafd-493c4cf70e84
# ╟─95c0b4f2-a01c-47b7-bb6e-64fbea5de550
# ╟─42fb46e3-3c37-483d-bafc-f4b9f4b11d1d
# ╟─b63d7fe7-845f-4f70-a5f8-81cfa5e1c018
# ╟─31eefb26-ac25-41f3-92ef-9a7600d6a093
# ╟─ab9084e1-57a2-4c6b-ad32-8fee8c142c43
# ╟─e6e1b8fd-05f5-4f86-a1d1-8a93f7a7854c
# ╠═37b61b8c-d983-4d67-938c-cd17b9749426
# ╠═5f7fb2f6-dbca-45ea-bb0a-1914260d0876
# ╟─18b29a1a-4787-11ec-25e3-5f29ebd21430
# ╟─16dca67c-f280-4a6f-bd79-308cf63dabf6
