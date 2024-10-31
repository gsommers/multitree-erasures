#!/usr/licensed/julia/1.10/bin/julia
module KnillPrep

using Base.Threads
using QuantumClifford
using StabilizerTree
using SteanePrep
using Utilities # sample_sites
using CircuitOps # trace_layer!
using TableauOps

export prep_measurements_reset, prep_measurements, prep_state_knill, get_bell_idxs, prep_knill_reset

function get_bell_idxs(t,tmax)
    gate_idxs = get_gate_idxs(t+1,tmax+1)
    cnot_idxs = vcat(gate_idxs[1:2:end], [[gate_idxs[j+1][1],gate_idxs[j][2]] for j=1:2:length(gate_idxs)])
    notc_idxs = vcat(gate_idxs[2:2:end], [[gate_idxs[j][1],gate_idxs[j+1][2]] for j=1:2:length(gate_idxs)])
    [cnot_idxs, notc_idxs]
end

# for the version where I reset qubits
# can be fixed n or not, with (possibly different) rate for each type of error
function get_erasure_sites_bell(rates, tmax)
    if length(rates)==1
        rates = fill(rates[1], 4)
    end
    input_sites = get_erasure_sites(tmax + 1, 2^tmax, rates[1]; fixed_n = false)

    # now I just need to fix up the first layer - where everything's an input
    input_sites = vcat([vcat(input_sites[1], input_sites[2] .+ 2^tmax)], input_sites[3:end])

    encoding_sites = get_erasure_sites(tmax, 2^(tmax+1), rates[2]; fixed_n = false)
    check_sites = get_erasure_sites(2*(tmax-1), 2^(tmax+1), rates[3]; fixed_n = false)
    check_sites2 = get_erasure_sites(tmax-1, 2^tmax, rates[3]; fixed_n = false)
    measure_sites = get_erasure_sites(tmax -1, 2^tmax, rates[4]; fixed_n = false)
    
    erasure_sites = [input_sites, encoding_sites, vcat(check_sites, check_sites2), measure_sites]
end

function prep_knill_reset(gate, tmax, rates, to_measure; input = :Z, erasure_sites = nothing)
    s, entropies, erasure_sites = prep_knill_reset(gate, tmax, rates; input = input, erasure_sites = erasure_sites)
    entropies[end] = system_entropy_knill!(copy(s), to_measure)
    s, entropies, erasure_sites
end

# version where I reset/feed back in measured qubits to keep a constant number
function prep_knill_reset(gate, tmax, rates; input = :Z, erasure_sites = nothing)
    s = tensor_pow(MixedDestabilizer(bell(1)⊗one(Stabilizer,2,basis=input)),2^(tmax-1))
    entropies = zeros(Int, tmax+1)

    if isnothing(erasure_sites)
        erasure_sites = get_erasure_sites_bell(rates, tmax)
    end
    big_gate = tensor_pow(gate, 2^(tmax))
    s_reset = one(Stabilizer, 2^tmax, basis=input)
    input_sites = 1:2^(tmax+1)

    # up to depth tmax...
    for t=1:tmax

    	gate_idxs = get_gate_idxs(t+1,tmax+1)
    	# Step 0: erasure errors on input legs
	trace_layer!(s, input_sites, erasure_sites[1][t])
	
        # Step 1: encoding gates
	apply!(s, big_gate, vcat(gate_idxs...))
	trace_layer!(s, 1:2^(tmax+1), erasure_sites[2][t])

	t==tmax && break
	# Step 2: CNOTs with ancillas
	bell_idxs = get_bell_idxs(t+1,tmax)
	for i=1:2
	    transversal_cnot!(s, bell_idxs[i], i; had = false)
	    trace_layer!(s, vcat(bell_idxs[i]...), erasure_sites[3][2*(t-1)+i])
	end
	
	# Step 3: Bell measurement of ancillas
	bell_pairs = [[bell_idxs[1][j][2],bell_idxs[2][j][2]] for j=1:length(bell_idxs[1])÷2]
	input_sites = vcat(bell_pairs...)
	# this way I'll measure half in the correct basis
	transversal_cnot!(s, bell_pairs, 2; had = true)
	trace_layer!(s, input_sites, erasure_sites[3][2*(tmax-1)+t])
	trace_layer!(s, input_sites, erasure_sites[4][t])
	for site in input_sites
	    projectZ!(s, site; keep_result = false)
	end
	entropies[t] = nqubits(s) - s.rank
	# reset inputs
	reset_qubits!(s, s_reset, input_sites)
    end
    entropies[tmax] = nqubits(s) - s.rank
    s, entropies, erasure_sites
end

function get_erasure_blocks(rates, tmax)
    encoding_sites = [[sample_sites(rates[2] * 2^(t+1), 2^(t+1); fixed_n = false) for i=1:3^(tmax-t)] for t=1:tmax]

    check_sites = [[sample_sites(rates[3] * 2^t, 2^t; fixed_n = false) for i=1:4*3^(tmax-t)] for t=2:tmax]

    measure_sites = [[sample_sites(rates[4] * 2^t, 2^t; fixed_n = false) for i=1:4*3^(tmax-t)] for t=2:tmax]

    input_sites = vcat([[sample_sites(rates[1] * 4, 4; fixed_n = false) for i=1:3^(tmax-1)]], [[sample_sites(rates[1] * 2^(t+1), 2^(t+1); fixed_n = false) for i=1:3^(tmax-t-1)] for t=1:tmax-1])

    [input_sites, encoding_sites, check_sites, measure_sites]
end

function make_bell(; input=:Z)
    s = bell(1)⊗one(Stabilizer,2,basis=input)
    MixedDestabilizer(s[:,[1,3,2,4]])
end

function system_entropy_knill!(s, to_measure)
    @assert nqubits(s)==nqubits(to_measure)
    for p in to_measure
        project!(s, p; keep_result = false)
    end
    ent = nqubits(s) - s.rank
    @assert ent<=2
    ent
end

function pair_blocks_knill(stabs, erasure_sites, to_measure; input = :Z)

    # versions where I do no checks
    entropies = [zeros(Int,3), zeros(Int,12), zeros(Int,12)]
    new_stabs = Array{MixedDestabilizer}(undef, 3)

    min_ent = 2
    best_idx = [1,1]
    for i=1:3
        new_stabs[i] = single_nocheck(stabs[i]; input = input)
	entropies[1][i] = system_entropy_knill!(copy(new_stabs[i]), to_measure)
	if entropies[1][i]==0
	    return new_stabs[i], [0,0], [1,i]
	end
	if entropies[1][i]<= min_ent
	    min_ent = entropies[1][i]
	    best_idx = [1,i]
	end
    end
    
    # versions where I do one check
    idxss = [[1,[1,1]],[1,[1,2]],[1,[2,1]],[1,[2,2]],[2,[1,1]],[2,[1,2]],[2,[2,1]],[2,[2,2]],[3,[1,1]],[3,[1,2]],[3,[2,1]],[3,[2,2]]]
    for i=1:length(idxss)
        new_stab = single_onecheck(stabs, [erasure_sites[1][1], erasure_sites[2][1]], idxss[i]; input = input)
	entropies[2][i] = system_entropy_knill!(new_stab, to_measure)
	if entropies[2][i]==0
	    # rerun to get the actual entropy
	    new_stab = single_onecheck(stabs, [erasure_sites[1][3], erasure_sites[2][3]], idxss[i]; input = input)
	    return new_stab,[0, system_entropy_knill!(copy(new_stab), to_measure)], [2,i]
	end
	if entropies[2][i]<= min_ent
	    min_ent = entropies[2][i]
	    best_idx = [2,i]
	end
    end

    # version where I do both checks
    # just default to the last one, for shits and grins
    for i=1:length(idxss)
        new_stab = single_twocheck(stabs, [erasure_sites[1][1:2], erasure_sites[2][1:2]], idxss[i]; input = input)
	entropies[3][i] = system_entropy_knill!(new_stab, to_measure)
	if entropies[3][i]==0
	    # rerun to get actual entropy
	    new_stab = single_twocheck(stabs, [erasure_sites[1][3:4], erasure_sites[2][3:4]], idxss[i]; input = input)
	    return new_stab, [0, system_entropy_knill!(copy(new_stab), to_measure)], [3,i]
	end
	if entropies[3][i]<= min_ent
	    min_ent = entropies[3][i]
	    best_idx = [3,i]
	end
    end
    if best_idx[1]==1 # best was "no check"
        return new_stabs[best_idx[2]], [min_ent, min_ent], best_idx
    elseif best_idx[1]==2 # best was "one check"
        new_stab = single_onecheck(stabs, [erasure_sites[1][3], erasure_sites[2][3]], idxss[best_idx[2]]; input = input)
        return new_stab,[min_ent, system_entropy_knill!(copy(new_stab), to_measure)], best_idx
    else
        new_stab = single_twocheck(stabs, [erasure_sites[1][3:4], erasure_sites[2][3:4]], idxss[best_idx[2]]; input = input)
	return new_stab, [min_ent, system_entropy_knill!(copy(new_stab), to_measure)], best_idx
    end
end

function prep_state_knill(gate, tmax, rates, to_measure; input = :Z, erasure_sites = nothing)
    n_blocks = 3^(tmax-1)
    my_s = make_bell(; input=input)

    stabs = [deepcopy(my_s) for i=1:n_blocks]
    entropies = [zeros(Int,2,3^(tmax-t-1)) for t=1:tmax-1]
    wirings = [zeros(Int,2,3^(tmax-t-1)) for t=1:tmax-1]

    if isnothing(erasure_sites)
        erasure_sites = get_erasure_blocks(rates, tmax)
    end
    
    big_gate = tensor_pow(gate,2)
    input_sites = 1:4
    gate_idxs = [[1,2],[3,4]]
    for t=1:tmax
        # Step 0: erasure errors on input legs
	for i=1:length(stabs)
            trace_layer!(stabs[i], input_sites, erasure_sites[1][t][i])
	end

        # Step 1: encoding gates
	for i=1:length(stabs)
	    apply!(stabs[i], big_gate, vcat(gate_idxs...))
	    trace_layer!(stabs[i], 1:2^(t+1), erasure_sites[2][t][i])
	end

	t==tmax && break
	big_gate = tensor_pow(big_gate,2)

	# Step 2: CNOTs with ancillas
	gate_idxs = get_gate_idxs(t+1,t+1)
	gate_idxs = vcat(gate_idxs, [idx .+ 2^(t+1) for idx in gate_idxs])
	input_sites = [idx[2] for idx in gate_idxs]
	
	new_stabs = Array{MixedDestabilizer}(undef, length(stabs)÷3)
	for i=1:length(stabs)÷3
	    new_stabs[i], entropies[t][:,i],wirings[t][:,i] = pair_blocks_knill(stabs[3*i-2:3*i], [erasure_sites[3][t][4*i-3:4*i], erasure_sites[4][t][4*i-3:4*i]], to_measure[t]; input = input)
	end
	stabs = new_stabs
    end

    final_ent = system_entropy_knill!(copy(stabs[end]), to_measure[end])

    stabs[end], final_ent, entropies, wirings, erasure_sites
end

function single_pair(stabs, erasure_sites, idxs; remove::Bool = true, input=:Z)
    # Knill measurement between (idxs[1]) half and (idxs[2]) half
    s = stabs[1]⊗stabs[2]

    n = nqubits(stabs[1])
    sites = [[i+(idxs[2]+1)*(n÷2),i+(idxs[1]-1)*(n÷2)] for i=1:nqubits(stabs[1])÷2]
    transversal_cnot!(s, sites, 2; had = true)

    use_sites = vcat(sites...)
    # errors associated with check gates
    trace_layer!(s, use_sites, erasure_sites[1])

    # errors associated with measurement
    trace_layer!(s, use_sites, erasure_sites[2])

    # measure
    for i=use_sites
        projectZ!(s, i; keep_result = false)
    end

    if remove
        # now cut out the part that was just measured
    	s = remove_qubits!(s, use_sites)
    else
        # reorder
	ordered_sites = vcat([i+(4-idxs[2])*(n÷2) for i=1:n÷2],[i+(idxs[1]-1)*(n÷2) for i=1:n÷2],[i+(2-idxs[1])*(n÷2) for i=1:n÷2],[i+(idxs[2]+1)*(n÷2) for i=1:n÷2])
	s.tab = s.tab[:,ordered_sites]
	if input !=:Z
            reset_qubits!(s, one(Stabilizer,n, basis=input), [n÷2+1:n;3*n÷2+1:2*n;])
	end
    end
    
    s
end

function single_nocheck(stab; input=:Z)
    s = stab⊗MixedDestabilizer(one(Stabilizer, nqubits(stab), basis=input))
    n = nqubits(s)
    s.tab = s.tab[:,[1:n÷4;n÷2+1:3*n÷4;n÷4+1:n÷2;3*n÷4+1:n;]]
    s
end

function single_onecheck(stabs, erasure_sites,idxs; input=:Z)
    use_blocks = setdiff(1:3,idxs[1])
    s = single_pair(stabs[use_blocks], erasure_sites, idxs[2]; remove = false, input=input)
end

function single_twocheck(stabs, erasure_sites, idxs; input = :Z)
    other_blocks = setdiff(1:3, idxs[1])
    s = single_pair([stabs[idxs[1]], stabs[other_blocks[1]]], [erasure_sites[1][1], erasure_sites[2][1]], [1,idxs[2][1]]; remove = true)
    s = single_pair([s, stabs[other_blocks[2]]], [erasure_sites[1][2], erasure_sites[2][2]], [1,idxs[2][2]]; remove = false, input = input)
end

function prep_measurements(gate, tmax; init_pauli = P"Z")
    to_measure = [s[:,get_system_idxs(t)] for (t,s)=enumerate(track_substabilizers(gate, tmax,init_pauli=init_pauli))]

    for t=1:tmax-1
        id_string = PauliOperator(zeros(Bool, 2*nqubits(to_measure[t])))
	to_measure[t] =tensor_pow(Stabilizer([p⊗id_string for p in to_measure[t]]),2)
    end
    to_measure[tmax] = tensor_pow(to_measure[tmax],2)

    to_measure
end

function prep_measurements_reset(gate, tmax; init_pauli = P"Z")
    s = Stabilizer(track_stabilizer_generators(gate, tmax, init_pauli = init_pauli))
    tensor_pow(s,2)[:,get_system_idxs(tmax+1)]
end

end # module