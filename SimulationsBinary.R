library(devtools)
library(latex2exp)
library(ggplot2)
library(scales)
library(reshape2)

set.seed('1234')

do_plot=FALSE
demo = FALSE

## DEMO ----
# just a short demo to show how the implemented classes work
if(demo){

  # let's generate some example stimuli input
  u = sample(0:1, 20, replace=T, prob=c(.5, .5))
  u = c(u, sample(0:1, 20, replace=T, prob=c(.2, .8)) )

  u = sample(0:1, 200, replace=T, prob=c(.2, .8))
  u = c(u, sample(0:1, 50, replace=T, prob=c(.8, .2)) )

  # now we can create a class with that input
  mathys = HGF_binary(u=u,parameters = list("kappa"=1.3, "omega"=-2.2, "theta"=0.5))
  any(mathys@moments$sigma < 0, na.rm=T)

  # the init function automatically calls mathys = fit(mathys, method="VB") if u is supplied
  plot(mathys)

  # the simulations have to be run explicitly, as they are time-consuming
  mathys = fit(mathys, method="RS")
  # this plots the empirical kernel estimation of distribution vs VB estimation
  plot.distributions(mathys, timestamps = c(10,50,80, 120, 150), levels = c(2,3))


  # the same goes for the continuous version
  u = c( rnorm(100, sd = 1) + 1:100/10, rnorm(100, sd=0.2)+10 )

  aaa = HGF_continuous(u=uu)
  plot(aaa)

  aaa = fit(aaa, method="BF")
  plot.distributions(aaa, timestamps = c(1,10,40,80), levels = c(1,2,3), logscale=T)

}

## SETTING UP STATISTICS FUNCTIONS ----


skewness = function(x){ mean( ( (x - mean(x))/sd(x) )**3 ) }

# estimate first three moments to save for later
compress_to_means = function(table){ apply(table, 2, mean) }
compress_to_vars = function(table){ apply(table, 2, var) }
compress_to_skewness = function(table){ apply(table, 2, skewness) }

# compress HGF instance to just the statistics from stochastic simulations
# (all the samples weight too much to store)
compressHGFInstance = function(instance){

  cnames = colnames(instance@simulations[[1]])

  means = sapply(instance@simulations, compress_to_means)
  means = t(means)
  colnames(means) = cnames

  vars = sapply(instance@simulations, compress_to_vars)
  vars = t(vars)
  colnames(vars) = cnames

  skew = sapply(instance@simulations, compress_to_skewness)
  skew = t(skew)
  colnames(skew) = cnames

  instance@simulations = list()
  instance@simulations[[1]] = means
  instance@simulations[[2]] = vars
  instance@simulations[[3]] = skew

  return(instance)
}

# accuracy statistic (whether the two learners agree)
get_dec_acc = function(instance){
  vb_dec = s(instance@moments[[1]][,2]) > 0.5
  mc_dec = instance@sim_preds_ > 0.5
  acc = mean(vb_dec == mc_dec)
  return(acc)
}

# smoothed relative error
maximal_sre = function(compressed_instance, levels=NA, moments=NA){
  if( any(is.na(levels)) ){
    levels = 2:compressed_instance@no_of_levels
    # since there is no noise, there is no need to check the first level
  }
  if( any(is.na(moments)) ){
    moments = 1:compressed_instance@no_of_moments
  }
  max_diff = 0
  for( level in levels ){
    for( moment in moments ){
      vb_moment = compressed_instance@moments[[moment]][,level]
      if( all(is.na(vb_moment)) ) next;
      mc_moment = compressed_instance@simulations[[moment]][,level]
      norm_dist = abs(mc_moment - vb_moment)/(abs(mc_moment) + 0.01)
      max_diff = max( max_diff, max(norm_dist, na.rm = T) )
    }
  }
 return(max(norm_dist))
}
average_sre = function(compressed_instance, levels=NA, moments=NA){
  if( any(is.na(levels)) ){
    levels = 2:compressed_instance@no_of_levels
    # since there is no noise, the first level will always be 0
  }
  if( any(is.na(moments)) ){
    moments = 1:compressed_instance@no_of_moments
  }
  diffs = 0
  for( level in levels ){
    for( moment in moments ){
      vb_moment = compressed_instance@moments[[moment]][,level]
      if( all(is.na(vb_moment)) ) next;
      mc_moment = compressed_instance@simulations[[moment]][,level]
      norm_dist = abs(mc_moment - vb_moment)/(abs(mc_moment) + 0.01)
      diffs = c( diffs, norm_dist )
    }
  }
  return(mean(diffs, na.rm = T))
}

pulled_sre = function(compressed_instance, levels=NA, moments=NA){
  if( any(is.na(levels)) ){
    levels = 2:compressed_instance@no_of_levels
    # since there is no noise, the first level will always be 0
  }
  if( any(is.na(moments)) ){
    moments = 1:compressed_instance@no_of_moments
  }
  diffs = 0
  divide_by = 0
  for( level in levels ){
    for( moment in moments ){
      vb_moment = compressed_instance@moments[[moment]][,level]
      if( all(is.na(vb_moment)) ) next;
      mc_moment = compressed_instance@simulations[[moment]][,level]
      norm_dist = abs(mc_moment - vb_moment)/(abs(mc_moment) + 0.01)
      diffs = diffs + norm_dist
      divide_by = divide_by + 1
    }
  }
  return(diffs/divide_by)
}

pulled_sre = function(compressed_instance, levels=NA, moments=NA){
  if( is.na(levels) ){
    levels = 2:compressed_instance@no_of_levels
    # since there is no noise, the first level will always be 0
  }
  if( is.na(moments) ){
    moments = 1:compressed_instance@no_of_moments
  }
  diffs = 0
  divide_by = 0
  for( level in levels ){
    for( moment in moments ){
      vb_moment = compressed_instance@moments[[moment]][,level]
      if( all(is.na(vb_moment)) ) next;
      mc_moment = compressed_instance@simulations[[moment]][,level]
      norm_dist = abs(mc_moment - vb_moment)/(abs(mc_moment) + 0.01)
      diffs = diffs + norm_dist
      divide_by = divide_by + 1
    }
  }
  return(diffs/divide_by)
}

## DEFINING SOME PLOT FUNCTIONS ----

plot_both_distributions = function(learner,
                                   actual_x_1=NA,
                                   point_size=0.2,
                                   title=NA,
                                   add_list=NA){

  if(is.na(title)){
    title = "$70\\%$ confidence intervals of $q_i$ and $\\hat{q}_i$ distributions"
  }

  N = length(learner@u)

  lvl1_df = data.frame(expand.grid(t=1:N,
                                   method=c("VB", "RS"),
                                   pred=NA))

  lvl1_df[ lvl1_df$method == "VB", "pred" ] =
    s( learner@moments$mu[,2] )
  lvl1_df[ lvl1_df$method == "RS", "pred" ] =
    sapply(learner@simulations, function(x) mean(s(x$x2)))

  lvl1 = ggplot() +
    geom_line(aes(x=t,
                  y=pred,
                  group=method,
                  colour=method), data=lvl1_df) +
    geom_point(aes(x=2:N,
                   y=learner@moments$mu[-1,1]),
               size = point_size) +
    ylab(TeX("$x_1$"))

  if( !all(is.na(actual_x_1)) ){
    lvl1 = lvl1 + geom_line(aes(x=2:N, y=actual_x_1))
  }

  # LEVEL 2
  lvl2_df = data.frame(expand.grid(t=1:N,
                                   method=c("VB", "RS"),
                                   mu=NA,
                                   lCI=NA,
                                   uCI=NA))
  lvl2_df[lvl2_df$method == "VB", "mu"] = learner@moments$mu[,2]
  lvl2_df[lvl2_df$method == "VB", "lCI"] = learner@moments$mu[,2] - sqrt(learner@moments$sigma[,2])*1
  lvl2_df[lvl2_df$method == "VB", "uCI"] = learner@moments$mu[,2] + sqrt(learner@moments$sigma[,2])*1

  lvl2_df[lvl2_df$method == "RS", c("lCI", "mu","uCI")] = t(sapply(learner@simulations,
                                                    function(x) quantile(x$x2, c(0.16,0.5,0.84)) ))

  lvl2 = ggplot(lvl2_df) +
    geom_line(aes(x=t,
                  y=mu,
                  group=method,
                  colour=method)) +
    geom_ribbon(aes(x=t,
                    ymin=lCI,
                    ymax=uCI,
                    group=method,
                    fill=method),
                alpha=0.2) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()) +
    ylab(TeX("$x_2$"))

  # LEVEL 3

  lvl3_df = data.frame(expand.grid(t=1:N,
                                   method=c("VB", "RS"),
                                   mu=NA,
                                   lCI=NA,
                                   uCI=NA))
  lvl3_df[lvl3_df$method == "VB", "mu"] = learner@moments$mu[,3]
  lvl3_df[lvl3_df$method == "VB", "lCI"] = learner@moments$mu[,3] - sqrt(learner@moments$sigma[,3])*1
  lvl3_df[lvl3_df$method == "VB", "uCI"] = learner@moments$mu[,3] + sqrt(learner@moments$sigma[,3])*1

  lvl3_df[lvl3_df$method == "RS", c("lCI", "mu","uCI")] = t(sapply(learner@simulations,
                                                                   function(x) quantile(x$x3, c(0.16,0.5,0.84)) ))

  lvl3 = ggplot(lvl3_df) +
    geom_line(aes(x=t,
                  y=mu,
                  group=method,
                  colour=method)) +
    geom_ribbon(aes(x=t,
                    ymin=lCI,
                    ymax=uCI,
                    group=method,
                    fill=method),
                alpha=0.2)+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()) +
    ylab(TeX("$x_3$")) +
    ggtitle(TeX(title))

  if( ! all(is.na(add_list)) ){
    if(length(add_list) == 1 ){
      for( element in add_list[1] ){ lvl1 = lvl1 + element }
      for( element in add_list[[1]] ){ lvl2 = lvl2 + element }
      for( element in add_list[[1]] ){ lvl3 = lvl3 + element }
    }else{
      for( element in add_list[[1]] ){ lvl1 = lvl1 + element }
      for( element in add_list[[2]] ){ lvl2 = lvl2 + element }
      for( element in add_list[[3]] ){ lvl3 = lvl3 + element }
    }
  }

  plot_grid(lvl3, lvl2, lvl1, ncol=1, align = "v")


}



## INITIAL PLOTS ----

set.seed(1234)
u = sample(0:1, 20, replace=T, prob=c(0.5,0.5))

learner1 = HGF_binary(u=u)

learner1 = fit(learner1, method = "RS")
figure_00 = plot_both_distributions(learner1, point_size = 1)
figure_00

cowplot::save_plot("figure00.jpg", figure_00, base_height = 4, dpi=400)

figure_01 = plot.distributions(learner1, timestamps = c(5, 10, 15, 20), levels = c(2,3))
figure_01
cowplot::save_plot("figure01.jpg", figure_01, base_height = 4, dpi=400)


set.seed(1234)

u = sample(0:1, 100, replace = T, prob=c(0.5,0.5))
u = c(u, sample(0:1, 20, replace = T, prob=c(0.8,0.2)))
u = c(u, sample(0:1, 20, replace = T, prob=c(0.2,0.8)))
u = c(u, sample(0:1, 20, replace = T, prob=c(0.8,0.2)))
u = c(u, sample(0:1, 20, replace = T, prob=c(0.2,0.8)))
u = c(u, sample(0:1, 20, replace = T, prob=c(0.8,0.2)))
u = c(u, sample(0:1, 100, replace = T, prob=c(0.5,0.5)))
true_x1 = c(rep(0.5, 100),
            rep(0.8, 20),
            rep(0.2, 20),
            rep(0.8, 20),
            rep(0.2, 20),
            rep(0.8, 20),
            rep(0.5, 100))
u = 1 - u

learner2 = HGF_binary(u=u)

learner2 = fit(learner2, method = "RS")
figure_02 = plot_both_distributions(learner2, actual_x_1 = true_x1)
figure_02

cowplot::save_plot("figure02.jpg", figure_02, base_height = 4, dpi=400)


set.seed(420)

u = sample(0:1, 100, replace = T, prob=c(0.5,0.5))
u = c(u, sample(0:1, 20, replace = T, prob=c(0.8,0.2)))
u = c(u, sample(0:1, 20, replace = T, prob=c(0.2,0.8)))
u = c(u, sample(0:1, 20, replace = T, prob=c(0.8,0.2)))
u = c(u, sample(0:1, 20, replace = T, prob=c(0.2,0.8)))
u = c(u, sample(0:1, 20, replace = T, prob=c(0.8,0.2)))
u = c(u, sample(0:1, 100, replace = T, prob=c(0.5,0.5)))
true_x1 = c(rep(0.5, 100),
            rep(0.8, 20),
            rep(0.2, 20),
            rep(0.8, 20),
            rep(0.2, 20),
            rep(0.8, 20),
            rep(0.5, 100))
u = 1 - u

learner3 = HGF_binary(u=u)

learner3 = fit(learner3, method = "RS")
figure_03 = plot_both_distributions(learner3, actual_x_1 = true_x1)
figure_03

cowplot::save_plot("figure03.jpg", figure_03, base_height = 4, dpi=400)


## SIMULATION 1 - BINARY STIMULI - SETTING UP DATA FRAME ----


if( file.exists("simTable1") ){
  simTable1=readRDS("simTable1")
}else{
  seed = 1:500
  p = seq(0.5,0.9,0.1)

  kappaSpace = 1.4
  omegaSpace = -2.2
  thetaSpace = 0.5

  DecisionAccuracy = NA

  completed = FALSE
  time_passed = NA

  simTable1 = expand.grid(
                          p=p,
                          seed=seed,
                          kappa=kappaSpace,
                          omega=omegaSpace,
                          theta=thetaSpace,
                          DecisionAccuracy=DecisionAccuracy,
                          completed = completed,
                          time_passed = time_passed)
}

# simArchives will be used to store the compressed HGF instances
if(file.exists("simArchives1")){
  simArchives1 = readRDS("simArchives1")
}else{
  simArchives1 = list()
}

### SIMULATION 1 - BINARY STIMULI - SIMULATION ----

# set simulation parameters
{
  time_rho1 = 50
  time_rho2 = 50
  max_time = time_rho1 + time_rho2
}

# set computation parameters (timers mostly)
{
  time_to_pass = 60*60*9
  save_interval = 60

  break_time = 60*15
  work_time = 60*60*4
}

# the simulation loop
{
  iter_no = which.min(simTable1$completed)

  # initialise the time variables
  start_time = Sys.time()
  end_time = start_time + time_to_pass
  last_save_time = Sys.time()
  last_break_time = Sys.time()

  # set up progress bar
  pb = txtProgressBar(min=0, max= nrow(simTable1), initial = sum(simTable1$completed), style=3)

  while( (iter_no <= nrow(simTable1)) & (Sys.time() < end_time) ){

    iteration_start_time = Sys.time()

    # read the parameters from the simulation table
    pars = as.list(simTable1[iter_no, c("kappa","omega","theta")])

    # generate the input stimuli "u"
    set.seed(simTable1[iter_no, "seed"])
    p = simTable1[iter_no, "p"]
    u = sample(0:1, time_rho1, replace=T, prob=c(1-p, p))
    u = c(u, sample(0:1, time_rho2, replace=T, prob=c(p, 1-p)) )

    # calculate both the HGF and the Rejection Sampling posteriors
    set.seed(iter_no) # to ensure different seed every simulation, but also replicability
    instance = HGF_binary(u=u,
                          parameters = pars)
    instance = fit(instance, method="RS")

    # store data in the table
    simTable1[iter_no,"DecisionAccuracy"] = get_dec_acc(instance)

    # store simulation time data
    simTable1[iter_no,"time_passed"] = as.numeric(Sys.time() - iteration_start_time , units="secs")

    # store the compressed instance in archive list
    simArchives1[[iter_no]] = compressHGFInstance(instance)

    # save the calculations in regular intervals
    if(as.numeric(Sys.time() - last_save_time , units="secs") > save_interval){
      saveRDS(simTable1,"simTable1")
      saveRDS(simArchives1,"simArchives1")
      last_save_time=Sys.time()
    }

    # update iteration number
    simTable1$completed[iter_no] = TRUE
    if(iter_no == nrow(simTable1)){
      iter_no = iter_no + 1
    }else{
      iter_no = which.min(simTable1$completed)
    }

    # update progress bar
    setTxtProgressBar(pb, sum(simTable1$completed))

    # sometimes cool down (I do this because I fear my weak computer may break)
    if(Sys.time() - last_break_time > work_time){
      last_break_time = Sys.time()
      Sys.sleep(break_time)
    }
  }
  close(pb)
}

### SIMULATION 1 - BIANRY STIMULI - PLOTTING RESULTS ----
if(do_plot){

  # Hit probability for a single p parameter
  p=0.8

  included = which( (simTable1$completed) & (simTable1$p==p) )

  plot_data_2 = data.frame(time=1:(max_time+1))
  for( i in included ){
    instance = simArchives1[[i]]
    vb_dec = s(instance@moments[[1]][,2]) >= 0.5
    mc_dec = instance@sim_preds_ >= 0.5

    plot_data_2[,paste("X", i)] = (vb_dec == mc_dec)
  }
  plot_data_2 = plot_data_2[, !apply(is.na(plot_data_2), 2, any)]
  plot_data_2$final = apply(plot_data_2[,-1], 1, mean)

  # calculating the confidence intervals
  alpha = 0.05
  z = qnorm(1 - alpha/2)
  n = length(included)
  CI_half = z/(sqrt(n)) * sqrt(0.25)
  plot_data_2$lower_CI = plot_data_2$final - z/(sqrt(n))*sqrt(plot_data_2$final*(1-plot_data_2$final))
  plot_data_2$upper_CI = plot_data_2$final + z/(sqrt(n))*sqrt(plot_data_2$final*(1-plot_data_2$final))

  plot_data_2$time = plot_data_2$time - min(plot_data_2$time)

  figure_1 = ggplot(plot_data_2[-1,]) +
    geom_line(aes(x=time, y=final)) +
    geom_ribbon(aes(x=time, ymin=lower_CI, ymax=upper_CI), alpha=0.15) +
    ggtitle(TeX(
      paste("Probability of getting the same prediction with $u\\sim Bern(", p, ")$")
      )) +
    ylab("")
  figure_1
  ggsave("figure1.jpg", figure_1, width = 15 , height = 10, units = "cm")


  # Hit probability for all parameters
  included_ps = unique(simTable1$p[ simTable1$completed ])

  plot_data_3 = expand.grid(time=1:(max_time+1), p=included_ps, hit_probability = NA)
  for( p in included_ps ){
    included = which( (simTable1$completed) & (simTable1$p==p) )
    helping_df = data.frame(time=1:(max_time+1))
    for( i in included ){
      instance = simArchives1[[i]]
      vb_dec = s(instance@moments[[1]][,2]) >= 0.5
      mc_dec = instance@sim_preds_ >= 0.5

      helping_df[,paste("X", i)] = (vb_dec == mc_dec)
    }
    helping_df = helping_df[, !apply(is.na(helping_df), 2, any)]
    helping_df$final = apply(helping_df[,-1], 1, mean)

    plot_data_3[ plot_data_3$p == p , "hit_probability" ] = helping_df$final

  }

  figure_2 = ggplot(plot_data_3) +
    geom_tile(aes(x = time, y=p, fill=hit_probability)) +
    ggtitle("Probability of getting the same prediction") +
    scale_y_continuous(breaks=unique(plot_data_3$p)) +
    scale_fill_binned(breaks=seq(0.45,1,0.05) , type = "viridis", direction=-1)+
    ylab(TeX("$\\rho_1$"))
  figure_2
  ggsave("figure2.jpg", figure_2, width = 15, height = 10, units = "cm")




  # Difference between moments for all parameters
  moment = 2
  level = 2
  included_ps = unique(simTable1$p[ simTable1$completed ])
  # included_ps = included_ps[1:3]
  plots_for_cow = list()
for( level in c(2,3)){
  plots_for_cow[[level]] = list()
  for(moment in c(1,2)){
    plot_data_3 = expand.grid(time=1:(max_time+1), p=included_ps, SRE = NA)
    for( p in included_ps ){
      included = which( (simTable1$completed) & (simTable1$p==p) )

      # kick out the outlier
      included = included[ included != 1650 ]

      helping_df = data.frame(time=1:(max_time+1))
      for( i in included ){
        instance = simArchives1[[i]]

        helping_df[,paste("X", i)] = pulled_sre(instance, level, moment)
      }
      helping_df = helping_df[, !apply(is.na(helping_df), 2, any)]
      helping_df$final = apply(helping_df[,-1], 1, mean)

      plot_data_3[ plot_data_3$p == p , "SRE" ] = helping_df$final

    }

    plots_for_cow[[level]][[moment]] =
    ggplot(plot_data_3) +
      geom_tile(aes(x = time, y=p, fill=SRE)) +
      scale_y_continuous(breaks=unique(plot_data_3$p)) +
      scale_fill_stepsn(limits=c(0,33),
                        breaks=c(seq(0,3,0.5),10,33),
                        colours=viridis::turbo(9),
                        values = scales::rescale(c(seq(0,3,0.5),10,33))) +
      ylab(TeX("$\\rho_1$"))
    if( moment == 1){

      plots_for_cow[[level]][[moment]] =plots_for_cow[[level]][[moment]] +
        ggtitle(TeX(
          paste("Averaged SRE($\\mu_" , level, ", \\hat{\\mu}_" , level , ")$")
        ))
    }else{

      plots_for_cow[[level]][[moment]] =plots_for_cow[[level]][[moment]] +
        ggtitle(TeX(
          paste("Averaged SRE($\\sigma_" , level, ", \\hat{\\sigma}_" , level , ")$")
        ))
    }
  }
}
cow = cowplot::plot_grid(plots_for_cow[[2]][[1]],
                   plots_for_cow[[2]][[2]],
                   plots_for_cow[[3]][[1]],
                   plots_for_cow[[3]][[2]],
                   ncol=2)

cow
cowplot::save_plot("figure3.jpg", cow, base_height = 5 )






# Now let's find out what generates the yellow mess..

plot(simArchives1[[1650]])

# Let's try to deal with the outliers

simTable1$pulled_SRE = NA

for( i in 1:nrow(simTable1)){
  simTable1$pulled_SRE[i] = average_sre(simArchives1[[i]])
}

simArchives1

}






## SIMULATION 2 - BINARY STIMULI - SETTING UP DATA FRAME ----


if( file.exists("simTable2") ){
  simTable2=readRDS("simTable2")
}else{
  seed = 1:25
  p = 0.75

  kappaSpace = seq(0,2, 0.2)
  omegaSpace = seq(-10,0, 1)
  thetaSpace = seq(0, 4, 0.4)

  DecisionAccuracy = NA

  completed = FALSE
  time_passed = NA

  simTable2 = expand.grid(
    p=p,
    kappa=kappaSpace,
    omega=omegaSpace,
    theta=thetaSpace,
    seed=seed,
    DecisionAccuracy=DecisionAccuracy,
    completed = completed,
    time_passed = time_passed)
}

# simArchives will be used to store the compressed HGF instances
if(file.exists("simArchives2")){
  simArchives2 = readRDS("simArchives2")
}else{
  simArchives2 = list()
}

### SIMULATION 2 - BINARY STIMULI - SIMULATION ----

# set simulation parameters
{
  time_rho1 = 50
  time_rho2 = 50
  max_time = time_rho1 + time_rho2
}

# set computation parameters (timers mostly)
{
  time_to_pass = 60*60*24*3
  save_interval = 60

  break_time = 60*15
  work_time = 60*60*4
}

# the simulation loop
{
  iter_no = which.min(simTable2$completed)

  # initialise the time variables
  start_time = Sys.time()
  end_time = start_time + time_to_pass
  last_save_time = Sys.time()
  last_break_time = Sys.time()

  # set up progress bar
  pb = txtProgressBar(min=0, max= nrow(simTable2), initial = sum(simTable2$completed), style=3)

  while( (iter_no <= nrow(simTable2)) & (Sys.time() < end_time) ){

    iteration_start_time = Sys.time()

    # read the parameters from the simulation table
    pars = as.list(simTable2[iter_no, c("kappa","omega","theta")])

    # generate the input stimuli "u"
    set.seed(simTable2[iter_no, "seed"])
    p = simTable2[iter_no, "p"]
    u = sample(0:1, time_rho1, replace=T, prob=c(1-p, p))
    u = c(u, sample(0:1, time_rho2, replace=T, prob=c(p, 1-p)) )

    # calculate both the HGF and the Rejection Sampling posteriors
    set.seed(iter_no) # to ensure different seed every simulation, but also replicability
    instance = HGF_binary(u=u,
                          parameters = pars)
    instance = fit(instance, method="RS")

    # store data in the table
    simTable2[iter_no,"DecisionAccuracy"] = get_dec_acc(instance)

    # store simulation time data
    simTable2[iter_no,"time_passed"] = as.numeric(Sys.time() - iteration_start_time , units="secs")

    # store the compressed instance in archive list
    simArchives2[[iter_no]] = compressHGFInstance(instance)

    # save the calculations in regular intervals
    if(as.numeric(Sys.time() - last_save_time , units="secs") > save_interval){
      saveRDS(simTable2,"simTable2")
      saveRDS(simArchives2,"simArchives2")
      last_save_time=Sys.time()
    }

    # update iteration number
    simTable2$completed[iter_no] = TRUE
    if(iter_no == nrow(simTable2)){
      iter_no = iter_no + 1
    }else{
      iter_no = which.min(simTable2$completed)
    }

    # update progress bar
    setTxtProgressBar(pb, sum(simTable2$completed))

    # sometimes cool down (I do this because I fear my weak computer may break)
    if(Sys.time() - last_break_time > work_time){
      last_break_time = Sys.time()
      Sys.sleep(break_time)
    }
  }
  close(pb)
}

### SIMULATION 2 - BIANRY STIMULI - PLOTTING RESULTS ----


if( do_plot ){


  simTable2 = simTable2[ simTable2$completed, ]

  simTable2$sre_mu_2 = sapply(simArchives2, function(x)average_sre(x, levels=2, moments=1) )
  simTable2$sre_mu_3 = sapply(simArchives2, function(x)average_sre(x, levels=3, moments=1) )
  simTable2$sre_sg_2 = sapply(simArchives2, function(x)average_sre(x, levels=2, moments=2) )
  simTable2$sre_sg_3 = sapply(simArchives2, function(x)average_sre(x, levels=3, moments=2) )
  names(simTable2)[9:12]
  simTable2$sre = apply(simTable2[,9:12], 1 , mean)

  plot_data_4 = simTable2[simTable2$sre < 1000, ]
  plot_data_4 = aggregate(plot_data_4[,6:13],
                          by=plot_data_4[,2:4],
                          FUN=function(x)mean(x,na.rm = T))

  plot_data_4 = plot_data_4[plot_data_4$omega==-2, ]
  plot_data_4 = plot_data_4[plot_data_4$kappa!=0, ]
  plot_data_4 = plot_data_4[plot_data_4$theta!=0, ]

  ggplot(plot_data_4) +
    geom_tile(aes(x=theta, y=kappa, fill=sre)) +
    facet_grid(omega~.) +
    scale_fill_continuous(limits=c(0,3),
                          type="viridis")

  figure21 = ggplot(plot_data_4[plot_data_4$theta == 2, ]) +
    geom_tile(aes(x=kappa, y=omega, fill=DecisionAccuracy)) +
    geom_text(aes(x=kappa, y=omega, label=round(DecisionAccuracy,2)), colour="white") +
    scale_fill_continuous(limits=c(0.5,1),
                          type="viridis",
                          direction=1) +
    scale_y_continuous(breaks = unique(plot_data_4$omega)) +
    scale_x_continuous(breaks = unique(plot_data_4$kappa))
  figure21
  ggsave2("figure21.jpg",figure21, width=6, height =4.5, dpi=400)


  plot_data_5 = plot_data_4
  plot_data_5 = rbind(plot_data_5,plot_data_5,plot_data_5,plot_data_5)
  nr = nrow(plot_data_5)
  plot_data_5$level = c(rep(2,nr/2), rep(3,nr/2))
  plot_data_5$moment = c(rep(1,nr/4), rep(2,nr/4), rep(1,nr/4), rep(2,nr/4))

  plot_data_5[(plot_data_5$moment == 1) & (plot_data_5$level == 2) , "sre" ] =
    plot_data_4$sre_mu_2
  plot_data_5[(plot_data_5$moment == 2) & (plot_data_5$level == 2) , "sre" ] =
    plot_data_4$sre_sg_2
  plot_data_5[(plot_data_5$moment == 1) & (plot_data_5$level == 3) , "sre" ] =
    plot_data_4$sre_mu_3
  plot_data_5[(plot_data_5$moment == 2) & (plot_data_5$level == 3) , "sre" ] =
    plot_data_4$sre_sg_3


  figure22 = ggplot(plot_data_5[plot_data_5$theta == 0.8,]) +
    geom_tile(aes(x=kappa, y=omega, fill=sre)) +
    facet_grid(moment~level) +
    scale_fill_continuous(limits=c(0,3),
                          type="viridis")+
    scale_y_continuous(breaks = unique(plot_data_4$omega)) +
    scale_x_continuous(breaks = unique(plot_data_4$kappa))
  figure22

  ggsave("figure22.jpg", figure22,width =6, height = 6, units="in", dpi=400)


  # INDIVIDUAL SIMUALTIONS

  recalculate = function( iter_no , time_rho1){

    pars = as.list(simTable2[iter_no, c("kappa","omega","theta")])
    set.seed(simTable2[iter_no, "seed"])
    p = simTable2[iter_no, "p"]
    u = sample(0:1, time_rho1, replace=T, prob=c(1-p, p))
    u = c(u, sample(0:1, time_rho2, replace=T, prob=c(p, 1-p)) )

    # calculate both the HGF and the Rejection Sampling posteriors
    set.seed(iter_no) # to ensure different seed every simulation, but also replicability
    instance = HGF_binary(u=u,
                          parameters = pars)
    instance = fit(instance, method="RS")

  }

  iter_no = which((simTable2$omega == -4 ) &
               (simTable2$kappa == 0.2 ) &
               (simTable2$theta == 0.8))[1]

  instance = recalculate(iter_no)

  fig23 = plot_both_distributions(instance)
  fig23
  ggsave("figure23.jpg", fig23,width =6, height = 6, units="in", dpi=400)


  {
    j=2
  iter_no = which((simTable2$omega == -8 ) &
                    (simTable2$kappa == 1 ) &
                    (simTable2$theta == 0.8))[j]
  instance = recalculate(iter_no)
  fig241 = plot_both_distributions(instance,
                                   title = "$\\hat{q}$ and $q$ distributions for $\\omega = -8$",
                                   add_list = list( list(theme(legend.position = "none")) ))

  iter_no = which((simTable2$omega == -4 ) &
                    (simTable2$kappa == 1 ) &
                    (simTable2$theta == 0.8))[j]
  instance = recalculate(iter_no)
  fig242 = plot_both_distributions(instance,
                                   title = "$\\hat{q}$ and $q$ distributions for $\\omega = -4$",
                                   add_list = list( list(theme(legend.position = "none"),
                                                         ylab("")) ))

  iter_no = which((simTable2$omega == -1 ) &
                    (simTable2$kappa == 1 ) &
                    (simTable2$theta == 0.8))[j]
  instance = recalculate(iter_no)
  fig243 = plot_both_distributions(instance,
                                   title = "$\\hat{q}$ and $q$ distributions for $\\omega = -1$",
                                   add_list = list( list(theme(legend.position = "none"),
                                                         ylab("")) ))

  fig24 = plot_grid(fig241, fig242, fig243, nrow=1)
  }
  fig24

  ggsave("figure24.jpg", fig24,width =12, height = 6, units="in", dpi=400)



  }



### end ----


if(do_plot){

  # Hit probability for a single p parameter
  p=0.8

  included = which( (simTable2$completed) & (simTable2$p-p<0.01) )

  plot_data_2 = data.frame(time=1:(max_time+1))
  for( i in included ){
    instance = simArchives2[[i]]
    vb_dec = s(instance@moments[[1]][,2]) >= 0.5
    mc_dec = instance@sim_preds_ >= 0.5

    plot_data_2[,paste("X", i)] = (vb_dec == mc_dec)
  }
  plot_data_2 = plot_data_2[, !apply(is.na(plot_data_2), 2, any)]
  plot_data_2$final = apply(plot_data_2[,-1], 1, mean)

  ggplot(plot_data_2[-1,]) +
    geom_line(aes(x=time, y=final)) +
    ggtitle(TeX(
      paste("Probability of getting the same prediction with $u\\sim Bern(", p, ")$")
    )) +
    ylab("")


  # Hit probability for all parameters
  included_ps = unique(simTable2$p[ simTable2$completed ])

  plot_data_3 = expand.grid(time=1:(max_time+1), p=included_ps, hit_probability = NA)
  for( p in included_ps ){
    included = which( (simTable2$completed) & (simTable2$p==p) )
    helping_df = data.frame(time=1:(max_time+1))
    for( i in included ){
      instance = simArchives2[[i]]
      vb_dec = s(instance@moments[[1]][,2]) >= 0.5
      mc_dec = instance@sim_preds_ >= 0.5

      helping_df[,paste("X", i)] = (vb_dec == mc_dec)
    }
    helping_df = helping_df[, !apply(is.na(helping_df), 2, any)]
    helping_df$final = apply(helping_df[,-1], 1, mean)

    plot_data_3[ plot_data_3$p == p , "hit_probability" ] = helping_df$final

  }

  ggplot(plot_data_3) +
    geom_tile(aes(x = time, y=p, fill=hit_probability)) +
    ggtitle("Probability of getting the same prediction") +
    scale_y_continuous(breaks=unique(plot_data_3$p)) +
    scale_fill_binned(breaks=(10:4)/10 , type = "viridis", direction=-1)+
    ylab("p")

  # Difference between moments for all parameters
  moment = 1
  level = 2
  included_ps = unique(simTable2$p[ simTable2$completed ])
  # included_ps = included_ps[1:3]

  plot_data_3 = expand.grid(time=1:(max_time+1), p=included_ps, distance = NA)
  for( p in included_ps ){
    included = which( (simTable2$completed) & (simTable2$p==p) )
    helping_df = data.frame(time=1:(max_time+1))
    for( i in included ){
      instance = simArchives2[[i]]
      vb_stat = instance@moments[[moment]][,level]
      mc_stat = instance@simulations[[moment]][,level]

      helping_df[,paste("X", i)] = abs(vb_stat - mc_stat)/abs(mc_stat)
    }
    helping_df = helping_df[, !apply(is.na(helping_df), 2, any)]
    helping_df$final = apply(helping_df[,-1], 1, mean)

    plot_data_3[ plot_data_3$p == p , "distance" ] = helping_df$final

  }

  ggplot(plot_data_3) +
    geom_tile(aes(x = time, y=p, fill=distance)) +
    ggtitle("Distance to true moment") +
    scale_y_continuous(breaks=unique(plot_data_3$p)) +
    ylab("diff")

}




### end ----


