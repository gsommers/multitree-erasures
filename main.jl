using SteaneRuns
using QuantumClifford

# example function you might call
function main_knill2(t,p)
    had2 = tHadamard⊗tHadamard
    gate = had2 * tCNOT
    @time entropy_dat = knill_erasures_reset([gate], t, [1], [[0,p,p,p]]; init_pauli="Z", num_samples=1000)
end
