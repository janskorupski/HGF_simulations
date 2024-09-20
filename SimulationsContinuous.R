library(devtools)
library(tuneR)
source("../HGFR/R/")

load_all("../HGFR/")

do_plot=FALSE

### DEMO ----
# just a short demo to show how the implemented classes work
if(T){

  u = c( rnorm(100, sd = 1) + 1:100/10, rnorm(100, sd=0.2)+10 )
  u = noise(kind = "pink", duration = 200)@left

  u = u*15
  set.seed(1234)
  u = noise(kind = "pink", duration = 200)@left
  aaa = HGF_continuous(u=u,
                       parameters = list("kappa"=rep(1.3,3), "omega"=rep(-2.2, 3), "theta"=4, "alpha"=0.01))
  any(aaa@moments$sigma < 0)
  plot(aaa)

  aaa = fit(aaa, method="BF")
  plot.distributions(aaa, timestamps = c(1,10,40,50,80,150), levels = c(1,2,3), logscale=T)

}

### Additional plot functions ----


plot_both_distributions = function(learner,
                                   actual_x_1=NA,
                                   title=NA,
                                   add_list=NA){

  if(is.na(title)){
    title = "$70\\%$ confidence intervals of $q_i$ and $\\hat{q}_i$ distributions"
  }

  N = length(learner@u)

  # LEVEL 1

  lvl1_df = data.frame(expand.grid(t=1:N,
                                   method=c("VB", "RS","u"),
                                   mu=NA,
                                   lCI=NA,
                                   uCI=NA))
  lvl1_df[lvl1_df$method == "VB", "mu"] = learner@moments$mu[,1]
  lvl1_df[lvl1_df$method == "VB", "lCI"] = learner@moments$mu[,1] - sqrt(learner@moments$sigma[,1])*1
  lvl1_df[lvl1_df$method == "VB", "uCI"] = learner@moments$mu[,1] + sqrt(learner@moments$sigma[,1])*1

  lvl1_df[lvl1_df$method == "RS", c("lCI", "mu","uCI")] = t(sapply(learner@simulations,
                                                                   function(x) weighted.quantiles(x$x1, x$w, c(0.16,0.5,0.84), na.rm=T) ))

  lvl1_df[lvl1_df$method == "u", "mu"] = learner@u
  lvl1_df[lvl1_df$method == "u", "lCI"] = NA
  lvl1_df[lvl1_df$method == "u", "uCI"] = NA

  lvl1_df = lvl1_df[(lvl1_df$method != "u") | (lvl1_df$t != 1),]

  lvl1 = ggplot(lvl1_df) +
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
    ylab(TeX("$x_1$"))



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
                                                                   function(x) weighted.quantiles(x$x2, x$w, c(0.16,0.5,0.84), na.rm=T) ))

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
                                                                   function(x) weighted.quantiles(x$x3, x$w, c(0.16,0.5,0.84), na.rm=T) ))

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

recalculate = function(iter_no, seed = 1234){

  set.seed(seed = seed)

  pars = as.list(simTable3[iter_no, c("kappa","omega","theta","alpha")])
  pars$kappa = rep(pars$kappa, 3)
  pars$omega = rep(pars$omega, 3)

  u = simArchives3[[iter_no]]@u[-1]

  instance = HGF_continuous(u=u,
                            parameters = pars)
  instance = fit(instance, method="BF")
}

### SIMULATION 3 - CONTINUOUS STIMULI - STIMULI GENERATION FUNCTION ----

generate_stimuli = function( max_time = 100,
                             parameters = list("kappa"=rep(1.4,3), "omega"=rep(-2.2, 3), "theta"=0.5, "alpha"=1),
                             priors = list("mu" = c(0, 0, 0), "sigma" = c(1, 1, 1))
){
  no_of_levels = length(parameters$kappa)
  u = NA

  x = matrix(NA, nrow=no_of_levels, ncol = max_time)

  for( level in 1:no_of_levels ){
    x[level,1] = rnorm(1,
                       mean = priors$mu[level],
                       sd = sqrt(priors$sigma[level]))
  }

  u = c(u, rnorm(1, mean = x[1,1], sd = sqrt(parameters$alpha)))

  for( t in 2:max_time ){
    x[no_of_levels,t] = rnorm(1,
                              mean = x[no_of_levels,t-1]*0.9,
                              sd = sqrt(parameters$theta) )
    for( level in (no_of_levels-1):1 ){
      var_upper = x[level + 1,t]*parameters$kappa[level] + parameters$omega[level]
      var_upper = exp(var_upper)
      x[level,t] = rnorm(1,
                         mean = x[level,t-1],
                         sd = sqrt(var_upper))
    }
    u = c( u, rnorm(1 ,
                    mean = x[1, t],
                    sd = sqrt(parameters$alpha)))
  }
  return( list(u=u, x=x) )
}


# SIMULATION 3 - CONTINUOUS STIMULI - STATISTIC AND HELP FUNCTIONS ----

# biased estimator, but it should not matter much for large sample sizes
weighted.variance = function(x, w, na.rm=T, average=NA){
  if(is.na(average)){
    average = weighted.mean(x, w, na.rm=na.rm)
  }
  variance = sum( w*( (x))**2 )/sum(w) - average**2
  return(variance)
}

# biased estimator as well
weighted.skewness = function(x, w, na.rm=T, average=NA){
  if(is.na(average)){
    average = weighted.mean(x, w, na.rm=na.rm)
  }
  variance = weighted.variance(x, w, na.rm=na.rm)
  skewness = sum( w*( (x - average)/sqrt(variance) )**3 )/sum(w)
  return(skewness)
}

weighted.quantiles = function(x, w, qantiles, na.rm=T){

  if( na.rm ){
    to_keep = (!is.na(x)) & (!is.na(w))
    w = w[ to_keep ]
    x = x[ to_keep ]
  }
  w = w/sum(w)

  new_idx = order(x)
  w = w[ new_idx ]
  x = x[ new_idx ]

  sw = cumsum(w)

  res=c()
  for( q in qantiles ){
    res = c(res, x[ which.min( sw < q ) ] )
  }
  return(res)

}

compress_to_means = function(table){

  table = table[ apply( !is.na(table), 1, all) , ]
  table = apply(table, 2, weighted.mean, w=table$w)
  return(table)
}

compress_to_vars = function(table){
  table = table[ apply( !is.na(table), 1, all) , ]
  table = apply(table, 2, weighted.variance, w=table$w)
}

compress_to_skew = function(table){
  table = table[ apply( !is.na(table), 1, all) , ]
  table = apply(table, 2, weighted.skewness, w=table$w)
}

compressHGFInstance = function(instance){

  cnames = colnames(instance@simulations[[1]])

  means = sapply(instance@simulations, compress_to_means)
  means = t(means)
  colnames(means) = cnames

  vars = sapply(instance@simulations, compress_to_vars)
  vars = t(vars)
  colnames(vars) = cnames

  skews = sapply(instance@simulations, compress_to_skew)
  skews = t(skews)
  colnames(skews) = cnames

  instance@simulations = list()
  instance@simulations$mu = means
  instance@simulations$sigma = vars
  instance@simulations$skew = skews

  return(instance)
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
      norm_dist = abs(mc_moment - vb_moment)/(abs(mc_moment) + 0.1)
      diffs = diffs + norm_dist
      divide_by = divide_by + 1
    }
  }
  return(diffs/divide_by)
}


### SIMULATION 3 - CONTINUOUS STIMULI - SETTING UP DATA FRAME ----

if( file.exists("simTable3") ){
  simTable3=readRDS("simTable3")
}else{
  seed=1:25

  kappaSpace = seq(0.4,2, 0.4)
  omegaSpace = seq(-10,0, 2)
  thetaSpace = seq(0.5,4,0.5)
  alphaSpace = seq(0.1,0.5,0.1)

  completed = FALSE
  time_passed = NA

  simTable3 = expand.grid(kappa=kappaSpace,
                          omega=omegaSpace,
                          theta=thetaSpace,
                          alpha=alphaSpace,
                          seed=seed,
                          completed = completed,
                          time_passed=time_passed)
}

if(file.exists("simArchives3")){
  simArchives3 = readRDS("simArchives3")
}else{
  simArchives3 = list()
}

### SIMULATIONS BINARY STIMULI - SIMULATION 1 ----

# set computation parameters (timers)
{
  time_to_pass = 60*60*20
  save_interval = 60*15

  break_time = 60*15
  work_time = 60*60*4
}

# initialise helping variables
{
  iter_no = which.min(simTable3$completed)

  start_time = Sys.time()
  end_time = start_time + time_to_pass
  last_save_time = Sys.time()
  last_break_time = Sys.time()

  pb = txtProgressBar(min=0, max= nrow(simTable3), initial = sum(simTable3$completed), style=3)
  while( (iter_no <= nrow(simTable3)) & (Sys.time() < end_time) ){

    iteration_start_time = Sys.time()

    pars = as.list(simTable3[iter_no, c("kappa","omega","theta","alpha")])
    pars$kappa = rep(pars$kappa, 3)
    pars$omega = rep(pars$omega, 3)

    u = noise(kind = "pink", duration = 100)@left

    instance = HGF_continuous(u=u,
                          parameters = pars)
    instance = fit(instance, method="BF")

    # store the simulations in archive list
    simArchives3[[iter_no]] = compressHGFInstance(instance)

    # save the calculations
    if(as.numeric(Sys.time() - last_save_time , units="secs") > save_interval){
      saveRDS(simTable3,"simTable3")
      saveRDS(simArchives3,"simArchives3")
      last_save_time=Sys.time()
    }

    iteration_end_time = Sys.time()
    # update iteration number
    simTable3$completed[iter_no] = TRUE
    simTable3$time_passed[iter_no] = iteration_end_time - iteration_start_time

    iter_no = which.min(simTable3$completed)

    setTxtProgressBar(pb, sum(simTable3$completed))

    if(Sys.time() - last_break_time > work_time){
      last_break_time = Sys.time()
      Sys.sleep(break_time)
    }
  }
  close(pb)
}

### PLOTTING RESULTS BINARY STIMULI - SIMULATION 1 ----
if(do_plot){

  # First let's calculate a bunch of interesting statistics

  simTable3 = simTable3[simTable3$completed, ]
  simArchives3 = simArchives3[ 1:(nrow(simTable3)) ]

  simTable3$SRE = sapply(simArchives3, average_sre)

  for( level in 1:3){
    simTable3[,paste("mu", as.character(level))] =
      sapply(simArchives3, function(x) average_sre( x, levels = level, moments = 1 ) )

    simTable3[,paste("sg", as.character(level))] =
      sapply(simArchives3, function(x) average_sre( x, levels = level, moments = 2 ) )
  }

  # SRE vs two variables

  plot_data_5 = simTable3
  plot_data_5 = aggregate(plot_data_5$SRE,
                          by=plot_data_5[,1:4],
                          FUN=function(x)mean(x,na.rm = T))

  fig31 = ggplot(plot_data_5) +
    geom_tile(aes(x=kappa,y=omega,fill=x)) +
    facet_grid(theta~alpha) +
    scale_y_continuous(breaks = unique(plot_data_5$omega),
                       sec.axis = sec_axis(~ . , name = "theta", breaks = NULL, labels = NULL)) +
    scale_x_continuous(breaks = unique(plot_data_5$kappa),
                       sec.axis = sec_axis(~ . , name = "alpha", breaks = NULL, labels = NULL))

  fig31
  ggsave2("figure31.jpg", fig31, width = 8, height = 8, units = "in", dpi=400)


  # plot 2
  plot_data_5 = simTable3
  plot_data_5 = aggregate(plot_data_5$SRE,
                          by=plot_data_5[,1:2],
                          FUN=function(x)mean(x,na.rm = T))

  fig32 = ggplot(plot_data_5) +
    geom_tile(aes(x=kappa,y=omega,fill=x)) +
    scale_y_continuous(breaks = unique(plot_data_5$omega)) +
    scale_x_continuous(breaks = unique(plot_data_5$kappa))

  fig32
  ggsave2("figure32.jpg", fig32, width = 8, height = 8, units = "in", dpi=400)

  # plot 3
  plot_data_6 = simTable3
  idxs = (abs(plot_data_6$theta - 2) < 0.01 )&
    (abs(plot_data_6$alpha - 0.3) < 0.01)
  plot_data_6 = plot_data_6[idxs , ]
  plot_data_6 = aggregate(plot_data_6[,9:14],
                          by=plot_data_6[,1:2],
                          FUN=function(x)mean(x,na.rm = T))

  plot_data_6 = melt( plot_data_6, id = c("kappa","omega") )
  plot_data_6$moment = substr(plot_data_6$variable, 1,2)
  plot_data_6$moment = c(mu=1,sg=2)[plot_data_6$moment]
  plot_data_6$level = substr(plot_data_6$variable, 4,4)
  plot_data_6$level = as.numeric(plot_data_6$level)
  plot_data_6$level = factor(plot_data_6$level, levels=3:1)
  names(plot_data_6)[4] = "SRE"

  fig33 = ggplot(plot_data_6) +
    geom_tile(aes(x = kappa, y = omega, fill = SRE)) +
    facet_grid(level~moment) +
    scale_y_continuous(breaks = unique(plot_data_5$omega),
                       sec.axis = sec_axis(~ . , name = "level", breaks = NULL, labels = NULL)) +
    scale_x_continuous(breaks = unique(plot_data_5$kappa),
                       sec.axis = sec_axis(~ . , name = "moment", breaks = NULL, labels = NULL))

  fig33
  ggsave2("figure33.jpg", fig33, width = 5, height = 6, units = "in", dpi=400)



  jj = which(   abs(simTable3$kappa - 0.4 ) < 0.1 &
                abs(simTable3$omega -  0  ) < 0.1 &
                abs(simTable3$theta - 2   ) < 0.1 &
                abs(simTable3$alpha - 0.3 ) < 0.1 )
  ll = recalculate(jj[7])
  fig34 = plot_both_distributions(ll)
  fig34

  ggsave("fig43.jpg", fig34, width=6, height = 6, dpi=400)


  # average pulled sre vs time
  levels = c(1,2,3)
  moments = c(2)

  idxs = (simTable3$kappa - 1.2 < 0.1) &
    (simTable3$omega == -10) &
    (simTable3$theta==0.5) &
    (simTable3$completed)
  idxs = which(idxs)
  sres =   sapply(simArchives3[idxs], pulled_sre, levels = levels, moments=moments)
  sres = apply(sres, 1, mean)

  base::plot(sres, type='l')

  # maximum sre vs time
  levels = c(1,2,3)
  moments = c(1)

  idxs = (simTable3$kappa - 1.2 < 0.1) &
    (simTable3$omega == -10) &
    (simTable3$theta==0.5) &
    (simTable3$completed)
  idxs = (simTable3$completed)
  idxs = which(idxs)
  sre =   sapply(simArchives3[idxs], pulled_sre, levels = levels, moments=moments)
  sre = apply(sre, 1, max)

  base::plot(sre, type='l')
}


# end ----
