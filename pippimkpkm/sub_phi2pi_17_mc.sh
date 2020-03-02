# ./create_sub.py "/lustre/nyx/panda/kgoetzen/gluex/data/simulation/phi2pi_17/tree_phi2pi_genr8_17.root" DSelector_pippimkpkm phi2pi_17_mc 200000 
# Data directory: $GLUEXDATA= /lustre/nyx/panda/kgoetzen/gluex/data/ 

sbatch job_rootmacro.sh 'run_dsel.C+("DSelector_pippimkpkm","simulation/phi2pi_17","pippimkpkm__B4","tree_phi2pi_genr8_17.root","phi2pi_17_mc_000000_000",200000,0)'  # ev = 200000
sbatch job_rootmacro.sh 'run_dsel.C+("DSelector_pippimkpkm","simulation/phi2pi_17","pippimkpkm__B4","tree_phi2pi_genr8_17.root","phi2pi_17_mc_000000_001",200000,200000)'  # ev = 200000
sbatch job_rootmacro.sh 'run_dsel.C+("DSelector_pippimkpkm","simulation/phi2pi_17","pippimkpkm__B4","tree_phi2pi_genr8_17.root","phi2pi_17_mc_000000_002",200000,400000)'  # ev = 200000
sbatch job_rootmacro.sh 'run_dsel.C+("DSelector_pippimkpkm","simulation/phi2pi_17","pippimkpkm__B4","tree_phi2pi_genr8_17.root","phi2pi_17_mc_000000_003",200000,600000)'  # ev = 200000
sbatch job_rootmacro.sh 'run_dsel.C+("DSelector_pippimkpkm","simulation/phi2pi_17","pippimkpkm__B4","tree_phi2pi_genr8_17.root","phi2pi_17_mc_000000_004",200000,800000)'  # ev = 200000
sbatch job_rootmacro.sh 'run_dsel.C+("DSelector_pippimkpkm","simulation/phi2pi_17","pippimkpkm__B4","tree_phi2pi_genr8_17.root","phi2pi_17_mc_000000_005",200000,1000000)'  # ev = 200000
sbatch job_rootmacro.sh 'run_dsel.C+("DSelector_pippimkpkm","simulation/phi2pi_17","pippimkpkm__B4","tree_phi2pi_genr8_17.root","phi2pi_17_mc_000000_006",200000,1200000)'  # ev = 200000
sbatch job_rootmacro.sh 'run_dsel.C+("DSelector_pippimkpkm","simulation/phi2pi_17","pippimkpkm__B4","tree_phi2pi_genr8_17.root","phi2pi_17_mc_000000_007",200000,1400000)'  # ev = 200000
sbatch job_rootmacro.sh 'run_dsel.C+("DSelector_pippimkpkm","simulation/phi2pi_17","pippimkpkm__B4","tree_phi2pi_genr8_17.root","phi2pi_17_mc_000000_008",200000,1600000)'  # ev = 200000
sbatch job_rootmacro.sh 'run_dsel.C+("DSelector_pippimkpkm","simulation/phi2pi_17","pippimkpkm__B4","tree_phi2pi_genr8_17.root","phi2pi_17_mc_000000_009",200000,1800000)'  # ev = 186027 

# EVENTS : 1,986,027
# JOBS   : 10
