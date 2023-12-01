# ------------------------------------------------------------------------

       # -- AutoNOM (Automatic Nonstationary Oscillatory Modelling) -- #

# ------------------------------------------------------------------------

using StatsBase
using Distributions
using Plots
using Optim
using DSP
using RCall
using ProgressMeter
using LinearAlgebra
using DelimitedFiles


path_stationary_fun = pwd()*"/functions_stationary_model.jl"
path_non_stationary_fun = pwd()*"/functions_non_stationary_model.jl"

include(path_stationary_fun)
include(path_non_stationary_fun)

global n_obs = 1000
Nsim = 50

# ------------ Parameters MCMC for AutoNOM ----------

n_iter_MCMC = 20000 # number of MCMC iteration for AutoNOM
n_CP_max = 15 # maximum number of change-points
n_freq_max = 10 # maximum number of frequencies in each segment.

# -- Segment model:

σ_β = 10 # prior variance for β,  β ∼ Normal(0, σ_β * I)
ν0 = 1/100 # prior σ², InverseGamma(ν0/2, η0/2)
γ0 = 1/100 # prior σ², InverseGamma(ν0/2, η0/2)
λ_S = 1  # poisson prior parameter, i.e n of freq ∼ Poisson(λ_S)
δ_ω_mixing_brn = 0.2 # mixing probability random walk when relocation a frequency
#                        (after 300 iterations burn-in)
c_S = 0.4 # constant for birth/death prbability, c ∈ (0, 0.5) -- Equation (6)
ϕ_ω = 0.25 # birth step, frequency is sampled from Unif(0, ψ_ω)
ψ_ω = 10/n_obs # miminum distance between frequencies
σ_RW_ω_hyper = 20 # variance parameter for random walk when relocating a freq :
#                   ωᵖ ∼ Normal(ωᶜ, σ_ω), where σ_ω = 1/(σ_RW_ω_hyper*n)

# -- Change-Point model:

λ_NS = 1/10 # poisson prior parameter, i.e n of cp ∼ Poisson(λ_NS)
c_NS = 0.4 # constant for birth/death prbability, c ∈ (0, 0.5) -- Equation (6)
ψ_NS = 10 # minumum distance between change-points
σ_RW_s = 2.5 # variance parameter for random walk when relocating a change-point
δ_s_mixing = 0.2 # mixing probability for mixture proposal when relocating a change-point

# ------------------------------------------------------ #
# ----------------- Simulation Study  ------------------ #
# ------------------------------------------------------ #
Nsim = 20
collect_cp = zeros(Nsim, 3)
collect_s = zeros(n_obs, Nsim)
collect_time = zeros(Nsim, 1)

for sim_index in 1:Nsim
  fname = join(["data_sim/setting_n1000_ncp6/data",sim_index,".txt"])

  wdf = open(fname,"r")

  global data = Array{Float64}(undef, n_obs, 1)
  header = readline(wdf, keep=false)
  index = 1

  while eof(wdf) != true
    tmp = readline(wdf, keep=true)
    data[index] = parse(Float64, tmp) 
    index += 1
  end

  # ----------- MCMC Objects and Starting Values -------------
  global y = copy(data)
  s_start = [40]
  MCMC_objects = get_MCMC_objects(s_start)

  β_sample = MCMC_objects["β"]
  ω_sample = MCMC_objects["ω"]
  σ_sample = MCMC_objects["σ"]
  s_sample = MCMC_objects["s"]
  n_freq_sample = MCMC_objects["n_freq"]
  log_likelik_sample = MCMC_objects["log_likelik"]
  n_CP_sample = MCMC_objects["n_CP"]

  # ------------- AutoNOM - Automatic Nonstationary Oscillatory Modelling  --------
  tstart = time()
  @showprogress for t in 2:(n_iter_MCMC + 1)

  global y = copy(data)
  global n = length(data)

  # Option: for the first 300 iterations,
  #         sample frequencies just from  periodogram (no Random Walk)

  if ( t <= 300 )
    global δ_ω_mixing = 1
  else
    global δ_ω_mixing = δ_ω_mixing_brn
  end

  n_CP_current = n_CP_sample[t-1]

  if (n_CP_current == 0)

    MCMC = birth_move_non_stationary(t, n_CP_sample, β_sample,
                  ω_sample, σ_sample, s_sample, n_freq_sample, log_likelik_sample)
    n_CP_sample = MCMC["n_CP"]
    β_sample = MCMC["β"]
    ω_sample = MCMC["ω"]
    σ_sample = MCMC["σ"]
    s_sample = MCMC["s"]
    n_freq_sample = MCMC["n_freq"]
    log_likelik_sample = MCMC["log_likelik"]

  elseif (n_CP_current ==  n_CP_max)

    death_prob_n_CP_max = c_NS*min(1, pdf(Poisson(λ_NS), n_CP_max-1)/pdf(Poisson(λ_NS), n_CP_max))
    within_prob_n_CP_max = 1 - death_prob_n_CP_max

    U = rand()

    if (U <= death_prob_n_CP_max)

      MCMC = death_move_non_stationary(t, n_CP_sample, β_sample,
                  ω_sample, σ_sample, s_sample, n_freq_sample, log_likelik_sample)
      n_CP_sample = MCMC["n_CP"]
      β_sample = MCMC["β"]
      ω_sample = MCMC["ω"]
      σ_sample = MCMC["σ"]
      s_sample = MCMC["s"]
      n_freq_sample = MCMC["n_freq"]
      log_likelik_sample = MCMC["log_likelik"]

    else

      MCMC = within_move_non_stationary(t, n_CP_sample, β_sample,
                  ω_sample, σ_sample, s_sample, n_freq_sample, log_likelik_sample)
      n_CP_sample = MCMC["n_CP"]
      β_sample = MCMC["β"]
      ω_sample = MCMC["ω"]
      σ_sample = MCMC["σ"]
      s_sample = MCMC["s"]
      n_freq_sample = MCMC["n_freq"]
      log_likelik_sample = MCMC["log_likelik"]

    end

  else

    # Probabilities for different types of moves
    birth_prob = c_NS*min(1, pdf(Poisson(λ_NS), n_CP_current+1)/pdf(Poisson(λ_NS), n_CP_current))
    death_prob = c_NS*min(1, pdf(Poisson(λ_NS), n_CP_current-1)/pdf(Poisson(λ_NS), n_CP_current))
    within_prob = 1 - (birth_prob + death_prob)

    U = rand()

    if (U <= birth_prob)

      MCMC = birth_move_non_stationary(t, n_CP_sample, β_sample,
                  ω_sample, σ_sample, s_sample, n_freq_sample, log_likelik_sample)
      n_CP_sample = MCMC["n_CP"]
      β_sample = MCMC["β"]
      ω_sample = MCMC["ω"]
      σ_sample = MCMC["σ"]
      s_sample = MCMC["s"]
      n_freq_sample = MCMC["n_freq"]
      log_likelik_sample = MCMC["log_likelik"]

    elseif ((U > birth_prob) && (U <= (birth_prob + death_prob)))

      MCMC = death_move_non_stationary(t, n_CP_sample, β_sample,
                  ω_sample, σ_sample, s_sample, n_freq_sample, log_likelik_sample)
      n_CP_sample = MCMC["n_CP"]
      β_sample = MCMC["β"]
      ω_sample = MCMC["ω"]
      σ_sample = MCMC["σ"]
      s_sample = MCMC["s"]
      n_freq_sample = MCMC["n_freq"]
      log_likelik_sample = MCMC["log_likelik"]


    else

      MCMC = within_move_non_stationary(t, n_CP_sample, β_sample,
                  ω_sample, σ_sample, s_sample, n_freq_sample, log_likelik_sample)
      n_CP_sample = MCMC["n_CP"]
      β_sample = MCMC["β"]
      ω_sample = MCMC["ω"]
      σ_sample = MCMC["σ"]
      s_sample = MCMC["s"]
      n_freq_sample = MCMC["n_freq"]
      log_likelik_sample = MCMC["log_likelik"]

    end

  end

  if (t % 2000 == 0)

    println("Iteration: ", t)

    println("Locations: ", s_sample[1:(n_CP_sample[t]), t])

    for j in 1:(n_CP_sample[t] + 1)

      println("Segment :", j)
      println("m: ", n_freq_sample[j, t])
      println("ω: ", ω_sample[j, 1:n_freq_sample[j, t], t])
      println("β: ", β_sample[j, 1:(2*n_freq_sample[j, t]+2), t])
      println("σ: ", σ_sample[j, t])

    end

  end
  end

  tend = time()


  # ----------------------- Diagnostic Convergence  ----------------------
  burn_in_MCMC = 5000
  final_indexes_MCMC = burn_in_MCMC:n_iter_MCMC

  # Trace log likelihood
  log_likelik_sample_total = zeros(Float64, n_iter_MCMC + 1)
  for t in 1:(n_iter_MCMC+1)
      n_CP = n_CP_sample[t]
      log_likelik_sample_total[t] = sum(log_likelik_sample[1:(n_CP+1), t])
  end

  # -- Markov Chains after burn-in
  n_CP_final = n_CP_sample[final_indexes_MCMC]
  β_final = β_sample[:, :, final_indexes_MCMC]
  ω_final = ω_sample[:, :, final_indexes_MCMC]
  s_final = s_sample[:, final_indexes_MCMC]
  σ_final = σ_sample[:, final_indexes_MCMC]
  n_freq_final = n_freq_sample[:, final_indexes_MCMC]
  n_final = length(n_CP_final)

  # Estimate: n of change-points
  n_CP_est = mode(n_CP_final)
  index_CP_est = findall(n_CP_final .== n_CP_est)

  # Estimate: n of frequencies, conditioned on n_CP_est
  n_freq_est = zeros(Int64, n_CP_est+1)
  for j in 1:(n_CP_est+1)
      n_freq_est[j] = mode(n_freq_final[j, index_CP_est])
  end

  # Estimate: change-points location
  s_est = Array{Float64,1}(undef,n_CP_est)
  for j in 1:n_CP_est
    s_est[j] = mean(s_final[j, index_CP_est])
  end

  collect_cp[sim_index,1:length(s_est)] = s_est

  # Estimate: signal, by averaging over MCMC iterations
  signal_MCMC = zeros(Float64, n_final, n_obs)
  @showprogress for t in 1:n_final
    n_CP = n_CP_final[t]
    signal_MCMC[t, :] = get_estimated_signal(data, s_final[1:n_CP, t],
                                      β_final[:, :, t],
                                      ω_final[:, :, t],
                                      n_freq_final[:, t])
  end

  sig_MCMC = zeros(Float64, n_obs, 1)
  for t in 1:n_obs
    sig_MCMC[t,1] = mean(signal_MCMC[:,t]) 
  end
  collect_s[:,sim_index] = sig_MCMC

  collect_time[sim_index] = tend - tstart

end

collect_time

collect_cp[1,:]
plot(data)
plot(collect_s[:,1])

writedlm("res_cp_julia.txt",  collect_cp, ',')
writedlm("res_time_julia.txt",  collect_time, ',')
writedlm("res_signal_julia.txt",  collect_s, ',')
