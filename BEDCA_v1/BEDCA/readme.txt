

This is for Euclidean distance correction of 3D spinning disc microscope

Use the driver - MasterRunExpts with some calling examples for real data, i.e. calls RunMCMC_eucldist, loading from library ExptData_list, or for simulated data call RunMCMC_ECsimulateddata (this loads library SimData_List)

EuclDistCorrection_Driver_fn (add new experiments in ExptData_List)

MCMC runs are saved in
MCMCruns
To also save angles use options in EuclDistCorrection_Driver_fn

You can load an mcmcrun with

load_mcmc_runs (unpackmcmcdat can be useful)

To generate statistics (automatic with MCMC run) use
create_outputs

Figures are:
Figures Posteriors for mu, sdx, sdz
MCMC_EuclDistMargFigs Posteriors of MC variables
MCMC_EuclDistMargFigsConv  Convergence plots
