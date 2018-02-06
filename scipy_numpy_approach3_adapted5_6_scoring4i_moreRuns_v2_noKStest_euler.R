# deSolve example
# "Chapter 3 Solving Ordinary Differential Equations in R"


library( "data.table" )
library( "deSolve" )
library( "ggplot2" )
library( "lhs" )
library( "MASS" )
library( "plyr" )
library( "reshape" )
library( "reshape2" )
library( "triangle" )
library( "truncnorm" )
library( "reshape" )


# setwd( "Q:/TB/Intrahost_model_prototypes/fitting_test" )
# setwd( "C:/Users/stchang/Documents" )
# setwd( "C:/Users/stchang/Documents/intrahost_candidate_models_111315d/reduce_incidence2c5" )
# setwd( "C:/Users/stewchang/Documents/EMOD/intrahost_candidate_models_111315d_versions/reduce_incidence2c5 - copy 170512 - seed 1 - fewer time windows - smaller priors - more runs redo" )
# setwd( "C:/Users/stewchang/Documents/EMOD/intrahost_candidate_models_111315d_versions/reduce_incidence2c5 - copy 170524 - new treg" )
setwd( "/home/ubuntu/treg07" )

flynn_data <- read.csv( "Lin_Flynn_2014_main_Fig4a.csv" )
green_data <- read.csv( "Green_2010_CD4_regulatory_T_cells.csv" )
green_data$Time <- green_data$Time * 7

temp_seed <- 7
write( temp_seed, file="temp_seed.txt" )
set.seed( temp_seed )

error_filehandle <- "sim_errors"

#====================================================
# Model parameters
#====================================================

# a <- -8/3
# b <- -10
# c <- 28
# yini <- c( X = 1, Y = 1, Z = 1 )


#----------------------------------------------------
# LHS 
#----------------------------------------------------

num_runs                <- 50000 # 5000 # 50000
num_params              <- 8
LHS.mat                 <- randomLHS( n=num_runs, k=num_params )

accept_rate             <- 0.025

# LHS.mat[,1] <- qunif( LHS.mat[,1], min=0, max=2 )
# LHS.mat[,2] <- qtriangle( LHS.mat[,2], a=0, b=10, c=1 )
# growth_rate_log_lhs     <- qunif( LHS.mat[,1], min=-4, max=-1 )
# imm_source_rate_log_lhs <- qunif( LHS.mat[,2], min=-8, max=-3 )

# priors <- cbind( growth_rate_log       = qunif( LHS.mat[,1], min=-4, max= -0.25 ),
                 # cmi_kill_rate_log     = qunif( LHS.mat[,2], min=-3, max= 0 ),
                 # imm_source_rate_log   = qunif( LHS.mat[,3], min=-8, max=-3 ),
                 # nai2pro_diff_rate_log = qunif( LHS.mat[,4], min=-3, max= 0 ),
                 # nai2ant_diff_rate_log = qunif( LHS.mat[,5], min=-3, max= 0 ),
                 # nai_decay_rate_log    = qunif( LHS.mat[,6], min=-4, max=-2 ),
                 # imm_ant_eff_log       = qunif( LHS.mat[,7], min=-1, max= 1 ) )
                 # # pro_decay_rate_log    = qunif( LHS.mat[,7], min=-4, max=-2 ),
                 # # ant_decay_rate_log    = qunif( LHS.mat[,8], min=-4, max=-2 ) )

# priors <- cbind( growth_rate_log       = qunif( LHS.mat[,1], min=-3, max= -0.25 ),
                 # cmi_kill_rate_log     = qunif( LHS.mat[,2], min=-3, max= 0 ),
                 # imm_source_rate_log   = qunif( LHS.mat[,3], min=-7, max=-2 ),
                 # nai2pro_diff_rate_log = qunif( LHS.mat[,4], min=-3, max= 0 ),
                 # nai2ant_diff_rate_log = qunif( LHS.mat[,5], min=-3, max= 0 ),
                 # nai_decay_rate_log    = qunif( LHS.mat[,6], min=-4, max=-1 ),
                 # pro_decay_rate_log    = qunif( LHS.mat[,7], min=-4, max=-1 ),
                 # ant_decay_rate_log    = qunif( LHS.mat[,8], min=-4, max=-1 ) )

# priors <- cbind(
                 # growth_rate_log       = qunif( LHS.mat[,1], min=-3, max= -0.25 ),
                 # # growth_rate_log       = qunif( LHS.mat[,1], min=-4, max= -0.25 ),
                 # cmi_kill_rate_log     = qunif( LHS.mat[,2], min=-3, max= 0 ),
                 # # imm_source_rate_log   = qunif( LHS.mat[,3], min=-5, max=-2 ),
                 # imm_source_rate_log   = qunif( LHS.mat[,3], min=-6, max=-2 ),
                 # nai2pro_diff_rate_log = qunif( LHS.mat[,4], min=-3, max= 0 ),
                 # nai2ant_diff_rate_log = qunif( LHS.mat[,5], min=-3, max= 0 ),
                 # nai_decay_rate_log    = qunif( LHS.mat[,6], min=-4, max=-2 ),
                 # pro_decay_rate_log    = qunif( LHS.mat[,7], min=-4, max=-2 ),
                 # ant_decay_rate_log    = qunif( LHS.mat[,8], min=-4, max=-2 ) )

priors <- cbind(
                 growth_rate_log       = qunif( LHS.mat[,1], min=-4, max= -0.25 ),
                 # growth_rate_log       = qunif( LHS.mat[,1], min=-3, max= -0.50 ),
                 # growth_rate_log       = qunif( LHS.mat[,1], min=-4, max= -0.25 ),
                 cmi_kill_rate_log     = qunif( LHS.mat[,2], min=-4, max= 0 ),
                 # imm_source_rate_log   = qunif( LHS.mat[,3], min=-5, max=-2 ),
                 imm_source_rate_log   = qunif( LHS.mat[,3], min=-7, max=-2 ),
                 nai2pro_diff_rate_log = qunif( LHS.mat[,4], min=-3, max= 0 ),
                 nai2ant_diff_rate_log = qunif( LHS.mat[,5], min=-3, max= 0 ),
                 nai_decay_rate_log    = qunif( LHS.mat[,6], min=-4, max=-2 ),
                 pro_decay_rate_log    = qunif( LHS.mat[,7], min=-4, max=-2 ),
                 ant_decay_rate_log    = qunif( LHS.mat[,8], min=-4, max=-2 ) )

# priors <- cbind( growth_rate_log       = qunif( LHS.mat[,1], min=-4, max=-0.25 ),
                 # cmi_kill_rate_log     = qunif( LHS.mat[,2], min=-3, max= 0 ),
                 # imm_source_rate_log   = qunif( LHS.mat[,3], min=-7, max=-4 ),
                 # nai2pro_diff_rate_log = qunif( LHS.mat[,4], min=-3, max= 0 ),
                 # nai2ant_diff_rate_log = qunif( LHS.mat[,5], min=-3, max= 0 ),
                 # nai_decay_rate_log    = qunif( LHS.mat[,6], min=-4, max=-2 ),
                 # pro_decay_rate_log    = qunif( LHS.mat[,7], min=-4, max=-2 ),
                 # ant_decay_rate_log    = qunif( LHS.mat[,8], min=-4, max=-2 ) )
rows   <- lapply( 1:num_runs, function(i){ priors[i,]} )


#----------------------------------------------------
# Bacterial parameters
#----------------------------------------------------

qui2act_precmi_rate     <- 3.0e-1    # 1.0e-3 # 0.1-0.5? *** Pre-immunity, rate of Mtb activation (greater) (poorly known, but conceivably measurable: can be estimated from measurements of EHR response?)
act2qui_precmi_rate     <- 3.0e-1    # 1.0e-3 # 0.001? *** Pre-immunity, rate of Mtb de-activation (assumed ratio of 1% inactive state)
growth_rate             <- 3.0e-1    # 1.0e-1
growth_rate_log         <- log10( growth_rate )
mtb_max_gra             <- 1.0e5     # 1.0e5 # 1.0e7? Maximum number of Mtb in a granuloma (carrying capacity) (well known from macaque model)
cmi_kill_rate           <- 6.0e-2    # 6.0e-2 # 1.0e-2 # 2.0e-1
cmi_kill_rate_log       <- log10( cmi_kill_rate )


#----------------------------------------------------
# Host parameters
#----------------------------------------------------

imm_source_rate         <- 1.0e-5    # 1.0e-3
imm_source_rate_log     <- log10( imm_source_rate )
nai2pro_diff_rate       <- 6.0e-1    # 1.0e-1 # 3.0e-2
nai2pro_diff_rate_log   <- log10( nai2pro_diff_rate )
nai2ant_diff_rate       <- 2.0e-2    # making this bigger makes activeation faster
nai2ant_diff_rate_log   <- log10( nai2ant_diff_rate )
imm_ant_eff             <- 10 # 1.0e0     # 1.0e0
imm_ant_eff_log         <- log10( imm_ant_eff )
nai_decay_rate          <- 1.0e-3    # 1.0e-4 # 1e-3 to 1e-2 is legit range, this is important
nai_decay_rate_log      <- log10( nai_decay_rate )
pro_decay_rate          <- 2.0e-3    # 2.0e0 # this is also important: 0 or 1 makes a difference
pro_decay_rate_log      <- log10( pro_decay_rate )
ant_decay_rate          <- 1.0e-3
ant_decay_rate_log      <- log10( ant_decay_rate )
# immunosen_rate          <- 3.0e-4    # 1.0e-4 # 5.0e-4 # 3.0e-4 # 1.0e-3
immunosen_rate          <- 0.0e0    # 1.0e-4 # 5.0e-4 # 3.0e-4 # 1.0e-3

cd4_max_bal <- 5.0e3    # Heron 2012


#----------------------------------------------------
# Prior distributions
#----------------------------------------------------

# priors <- cbind( growth_rate_log       = qtruncnorm( LHS.mat[,1], a=growth_rate_log-1, b=growth_rate_log+1, mean=growth_rate_log, sd=0.5 ),
                 # cmi_kill_rate_log     = qtruncnorm( LHS.mat[,2], a=cmi_kill_rate_log-1, b=cmi_kill_rate_log+1, mean=cmi_kill_rate_log, sd=0.5 ),
                 # imm_source_rate_log   = qtruncnorm( LHS.mat[,3], a=imm_source_rate_log-1, b=imm_source_rate_log+1, mean=imm_source_rate_log, sd=0.5 ),
                 # nai2pro_diff_rate_log = qtruncnorm( LHS.mat[,4], a=nai2pro_diff_rate_log-1, b=0, mean=nai2pro_diff_rate_log, sd=0.5 ),
                 # nai2ant_diff_rate_log = qtruncnorm( LHS.mat[,5], a=nai2ant_diff_rate_log-1, b=0, mean=nai2ant_diff_rate_log, sd=0.5 ),
                 # nai_decay_rate_log    = qtruncnorm( LHS.mat[,6], a=nai_decay_rate_log-1, b=nai_decay_rate_log+1, mean=nai_decay_rate_log, sd=0.5 ),
                 # pro_decay_rate_log    = qtruncnorm( LHS.mat[,7], a=pro_decay_rate_log-1, b=pro_decay_rate_log+1, mean=pro_decay_rate_log, sd=0.5 ),
                 # ant_decay_rate_log    = qtruncnorm( LHS.mat[,8], a=ant_decay_rate_log-1, b=ant_decay_rate_log+1, mean=ant_decay_rate_log, sd=0.5 ) )
# rows   <- lapply( 1:num_runs, function(i){ priors[i,]} )


#----------------------------------------------------
# Model initial conditions
#----------------------------------------------------

### Bacterial initial conditions ###
mtb_qui_gra_init        <- 0
mtb_act_gra_init        <- 1                 # Initial condition of latent disease is 1 active Mtb bacillus

### Host initial conditions ###
imm_nai_init            <- 1                 # CD4 count
imm_pro_init            <- 0
imm_ant_init            <- 0


#----------------------------------------------------
# Model steps
#----------------------------------------------------

start_time <- 0
end_time   <- 700
intervals  <- 100 # 1000

timesteps <- seq( from=start_time, to=end_time, length.out=intervals )


#====================================================
# Model equations
#====================================================

model_eq <- function( timepoint, y_at_timepoint, parms ) {

    # def funct( state_vars, time_step ):
    with( as.list( c( y_at_timepoint, parms ) ), {
        
        # # state variables        
        # mtb_qui_gra_t = state_vars[0]
        # mtb_act_gra_t = state_vars[1]
        # imm_nai_t     = state_vars[2]
        # imm_pro_t     = state_vars[3]
        # imm_ant_t     = state_vars[4]
        
        # variable parameters
        # growth_rate_log, imm_source_rate_log, nai2pro_diff_rate_log, nai2ant_diff_rate_log, nai_decay_rate_log, pro_decay_rate_log, ant_decay_rate_log = parms
        
        # model equations
        
        mtb_qui_gra_delta  <- - ( qui2act_precmi_rate * mtb_qui_gra_t ) + ( act2qui_precmi_rate * mtb_act_gra_t )
        ###                  - activation (with effector-T regulation) + de-activation ###
        
        # mtb_act_gra_delta  <-   ( qui2act_precmi_rate * mtb_qui_gra_t ) - ( act2qui_precmi_rate * mtb_act_gra_t ) + 10 ** growth_rate_log * mtb_act_gra_t * ( 1 - ( mtb_qui_gra_t + mtb_act_gra_t ) / mtb_max_gra ) - ( ( 10 ** cmi_kill_rate_log ) * max( 0, ( imm_pro_t - imm_ant_eff * imm_ant_t ) ) * mtb_act_gra_t ) # - ( ipt_effect * mtb_act_gra_t )
        # mtb_act_gra_delta  <-   ( qui2act_precmi_rate * mtb_qui_gra_t ) - ( act2qui_precmi_rate * mtb_act_gra_t ) + 10 ** growth_rate_log * mtb_act_gra_t * ( 1 - ( mtb_qui_gra_t + mtb_act_gra_t ) / mtb_max_gra ) - ( ( 10 ** cmi_kill_rate_log ) * imm_pro_t * ( imm_pro_t / ( imm_nai_t + imm_ant_eff * imm_ant_t ) ) * mtb_act_gra_t ) # - ( ipt_effect * mtb_act_gra_t )
        mtb_act_gra_delta  <-   ( qui2act_precmi_rate * mtb_qui_gra_t ) - ( act2qui_precmi_rate * mtb_act_gra_t ) + ( 10 ** growth_rate_log ) * mtb_act_gra_t * ( 1 - ( mtb_qui_gra_t + mtb_act_gra_t ) / mtb_max_gra ) - ( 10 ** cmi_kill_rate_log ) * max( 0, ( imm_pro_t - imm_ant_t ) ) * mtb_act_gra_t # - ( ipt_effect * mtb_act_gra_t )
        ###                     activation (with effector-T regulation) - de-activation                           + replication                                                                                     - cmi-mediated killing ###
        
        # imm_nai_delta      <-   ( 10 ** imm_source_rate_log ) * ( ( 1 - immunosen_rate / 365. ) ** timepoint ) * mtb_qui_gra_t - ( 10 ** nai2pro_diff_rate_log * imm_nai_t ) - ( 10 ** nai2ant_diff_rate_log * imm_nai_t ) - ( 10 ** nai_decay_rate_log * imm_nai_t )
        # imm_nai_delta      <-   ( 10 ** imm_source_rate_log ) * ( ( 1 - immunosen_rate / 365. ) ** timepoint ) * mtb_qui_gra_t - ( 10 ** nai2pro_diff_rate_log * imm_nai_t ) - ( 10 ** nai2ant_diff_rate_log * imm_nai_t ) - ( 10 ** nai_decay_rate_log * imm_nai_t )
        # imm_nai_delta      <-   ( 10 ** imm_source_rate_log ) * ( ( mtb_qui_gra_t + mtb_act_gra_t ) / mtb_max_gra ) * ( 1 - ( imm_nai_t + imm_pro_t + imm_ant_t ) / cd4_max_bal ) - ( 10 ** nai2pro_diff_rate_log * imm_nai_t ) - ( 10 ** nai2ant_diff_rate_log * imm_nai_t ) - ( 10 ** nai_decay_rate_log * imm_nai_t )
        imm_nai_delta      <-   ( 10 ** imm_source_rate_log ) * mtb_act_gra_t * ( 1 - ( mtb_qui_gra_t + mtb_act_gra_t ) / mtb_max_gra ) * ( 1 - ( imm_nai_t + imm_pro_t + imm_ant_t ) / cd4_max_bal ) - ( 10 ** nai2pro_diff_rate_log * imm_nai_t ) - ( 10 ** nai2ant_diff_rate_log * imm_nai_t ) - ( 10 ** nai_decay_rate_log * imm_nai_t )
        ###                     T-cell replication                                                                                                                   - effector-T differentiation                  - regulatory-T differentiation                - decay ###
        
        imm_pro_delta      <-                                                                                                                                          ( 10 ** nai2pro_diff_rate_log * imm_nai_t )                                               - ( 10 ** pro_decay_rate_log * imm_pro_t )
        ###                                                                                                                                                            effector-T differentiation                                                                - decay ###
        
        imm_ant_delta      <-                                                                                                                                                                                        ( 10 ** nai2ant_diff_rate_log * imm_nai_t ) - ( 10 ** ant_decay_rate_log * imm_ant_t )
        ###                                                                                                                                                                                                          regulatory-T differentiation                - decay ###
        
        return( list( c( mtb_qui_gra_delta, mtb_act_gra_delta, imm_nai_delta, imm_pro_delta, imm_ant_delta ) ) )
        
    } )
    
}


#====================================================
# Model simulations
#====================================================

y0         <- c( mtb_qui_gra_t=mtb_qui_gra_init, mtb_act_gra_t=mtb_act_gra_init, imm_nai_t=imm_nai_init, imm_pro_t=imm_pro_init, imm_ant_t=imm_ant_init )

# init_rates <- c( growth_rate_log=growth_rate_log, cmi_kill_rate_log=cmi_kill_rate_log, imm_source_rate_log=imm_source_rate_log, nai2pro_diff_rate_log=nai2pro_diff_rate_log, nai2ant_diff_rate_log=nai2ant_diff_rate_log, nai_decay_rate_log=nai_decay_rate_log, pro_decay_rate_log=pro_decay_rate_log, ant_decay_rate_log=ant_decay_rate_log )

init_rates <- c(    growth_rate_log=growth_rate_log,
                    cmi_kill_rate_log=cmi_kill_rate_log,
                    imm_source_rate_log=imm_source_rate_log,
                    nai2pro_diff_rate_log=nai2pro_diff_rate_log,
                    nai2ant_diff_rate_log=nai2ant_diff_rate_log,
                    nai_decay_rate_log=nai_decay_rate_log,
                    pro_decay_rate_log=pro_decay_rate_log,
                    ant_decay_rate_log=ant_decay_rate_log )

# init_rates <- c( growth_rate_log=growth_rate_log, cmi_kill_rate_log=cmi_kill_rate_log, imm_source_rate_log=imm_source_rate_log, nai2pro_diff_rate_log=nai2pro_diff_rate_log, nai2ant_diff_rate_log=nai2ant_diff_rate_log, nai_decay_rate_log=nai_decay_rate_log, imm_ant_eff_log=imm_ant_eff_log )

# Sample parameters and use in solving ODE's
# sims       <- lapply( rows, function( lhs_params ){ init_rates[ names( lhs_params ) ] <- lhs_params
                                                    # print( init_rates )
                                                    # ode( y=y0, times=timesteps, func=model_eq, parms=init_rates ) } )

write( names( init_rates ), file=paste( error_filehandle, ".csv", sep="" ), ncolumns=length(init_rates), sep=",", append=F )
sims       <- llply( .data = rows, .fun = function( params ){
                                                                init_rates[ names( params[ 1:length(init_rates) ] ) ] <- params
                                                                print( init_rates )
                                                                # ode( y=y0, times=timesteps, func=model_eq, parms=init_rates, rtol=1e-15, method="rk4" )
                                                                # tryCatch( ode( y=y0, times=timesteps, func=model_eq, parms=init_rates, rtol=1e-15, maxsteps=20000, method="lsoda" ),
                                                                tryCatch( ode( y=y0, times=timesteps, func=model_eq, parms=init_rates, method="euler", hini=1 ),
                                                                ### example: ode( y=state, times=0:40, func=Lorenz, parms=parameters, method="euler", hini=0.01)
                                                                warning = function(w){ print(w); write( init_rates, file=paste( error_filehandle, ".csv", sep="" ), ncolumns=length(init_rates), sep=",", append=T ) } )
                                                            }, .progress= progress_text(char = '*') )


# # Try stiff solver?
# sims       <- llply( .data = rows, .fun = function( lhs_params ){ init_rates[ names( lhs_params ) ] <- lhs_params
                                                                  # # print( init_rates )
                                                                  # lsoda( y=y0, times=timesteps, func=model_eq, parms=init_rates )
                                                                # }, .progress= progress_text(char = '*') )

dead_sims <- which( sapply( sims, is.null ) )


#====================================================
# Model scoring
#====================================================

# flynn_time_divs  <- list( c( 0, 50 ), c( 50, 200 ), c( 200, 300 ), c( 300, 400 ), c( 400, 600 ) )
# flynn_time_divs  <- list( c( 0, 50 ), c( 50, 200 ), c( 200, 400 ), c( 400, 600 ) )
# flynn_time_divs  <- list( c( 0, 200 ), c( 200, 400 ), c( 400, 600 ) )
# flynn_time_divs  <- list( c( 0, 50 ), c( 50, 200 ), c( 200, 600 ) )
flynn_time_divs  <- list( c( 0, 120 ), c( 120, 300 ), c( 300, 600 ) )
# flynn_time_divs  <- list( c( 0, 115 ), c( 115, 600 ) )

# flynn_time_divs  <- list( c( 0, 50 ), c( 50, 150 ), c( 150, 350 ), c( 350, 600 ) )
# flynn_time_divs  <- list( c( 0, 50 ), c( 50, 300 ), c( 300, 600 ) )
# flynn_time_divs  <- list( c( 0, 100 ), c( 100, 300 ), c( 300, 600 ) )

# green_time_divs  <- list( c( 0, 50 ), c( 50, 100 ), c( 100, 150 ), c( 150, 170 ) )
green_time_divs  <- list( c( 0, 50 ), c( 50, 100 ), c( 100, 170 ) )
# green_time_divs  <- list( c( -1, 50 ), c( 50, 125 ), c( 125, 170 ) )

# flynn_time_divs  <- list( c( -1, 100 ), c( 100, 300 ), c( 300, 600 ) )
# green_time_divs  <- list( c( -1, 50 ), c( 50, 150 ), c( 150, 170 ) )

calc_stat_sum <- function( sims_df, time_divs, data_set, data_var ) {
    
    # stat_sum <- stat_sum_norm <- stat_sum_alt <- stat_sum_z <- 0
    stat_sum <- stat_sum_alt <- stat_sum_z <- 0
    # cfu_stat_sum  <- cfu_stat_sum_norm  <- cfu_stat_sum_alt  <- 0
    # treg_stat_sum <- treg_stat_sum_norm <- treg_stat_sum_alt <- 0
    
    for( j in seq_along( time_divs ) ) {
        
        # ---------- Scoring for CFU data ----------
        if( tolower( data_var ) == "cfu" ) {
            
            # ----- Stat #1: KS test statistic ----- #
            # ks_res        <- ks.test( x = log10( sims_df [ ( ( sims_df$time  > time_divs[[j]][1] ) & ( sims_df$time  < time_divs[[j]][2] ) ), ]$mtb_act_gra_t + 1 ),
                                      # y = log10( data_set[ ( ( data_set$Time > time_divs[[j]][1] ) & ( data_set$Time < time_divs[[j]][2] ) ), ]$CFU ),
                                      # alternative = "t", exact = FALSE )
            
            # # ----- Stat #2: KS test statistic normalized ----- #
            # ks_stat_norm  <- ks_res$statistic / sd( log10( data_set[ ( data_set$Time > time_divs[[j]][1] ) & ( data_set$Time < time_divs[[j]][2] ), ]$CFU ) )
            
            # ----- Stat #3: Z-score ----- #
            sims_sub_mu   <- mean( log10( sims_df [ ( ( sims_df$time  > time_divs[[j]][1] ) & ( sims_df$time  < time_divs[[j]][2] ) ), ]$mtb_act_gra_t + 1 ) )
            data_sub_mu   <- mean( log10( data_set[ ( ( data_set$Time > time_divs[[j]][1] ) & ( data_set$Time < time_divs[[j]][2] ) ), ]$CFU ) )
            data_sub_sd   <- sd  ( log10( data_set[ ( ( data_set$Time > time_divs[[j]][1] ) & ( data_set$Time < time_divs[[j]][2] ) ), ]$CFU ) )
            
            # alt_stat      <- abs( sims_sub_mu - data_sub_mu ) / ( data_sub_sd / sqrt( length( log10( data_set[ ( ( data_set$Time > time_divs[[j]][1] ) & ( data_set$Time < time_divs[[j]][2] ) ), ]$CFU ) ) ) )         # Z-score
            alt_stat      <- abs( sims_sub_mu - data_sub_mu ) / ( data_sub_sd )         # Z-score
            # print( alt_stat )
            
            # ----- Stat #4: Z-score 2 ----- #
            
            sim_subset     <- sims_df[ ( ( sims_df$time  > time_divs[[j]][1] ) & ( sims_df$time < time_divs[[j]][2] ) ), ]
            data_subset    <- data_set[ ( ( data_set$Time > time_divs[[j]][1] ) & ( data_set$Time < time_divs[[j]][2] ) ), ]
            data_subset_sd <- sd( log10( data_subset$CFU ) )
            
            z_score <- 0
            for( k in 1:nrow( data_subset ) ) {
                # print( data_subset$Time[k] )
                pos_time_match <- which.min( abs( data_subset$Time[k] - sim_subset$time ) )
                z_score <- z_score + abs( log10( data_subset$CFU[k] ) - log10( sim_subset$mtb_act_gra_t[pos_time_match] + 1 ) ) / data_subset_sd        # define a score as the average absolute value of z-score
            }
            z_score <- z_score / nrow( data_subset )
            
        # ---------- Scoring for Treg data ----------
        } else if( tolower( data_var ) == "tregpop" ) {
            
            sims_sub      <- sims_df[ ( sims_df$time  > time_divs[[j]][1] ) & ( sims_df$time  < time_divs[[j]][2] ), ]
            sims_sub_perc <- 100 * sims_sub$imm_ant_t / ( sims_sub$imm_nai_t + sims_sub$imm_pro_t + sims_sub$imm_nai_t )
            # print( c( min( sims_sub_perc ), max( sims_sub_perc ) ) )
            # print( sims_sub_perc )
            
            sims_sub_perc_mu <- mean( sims_sub_perc )
            # sims_sub_perc_sd <- sd  ( sims_sub_perc$Mean )
            
            green_data_N  <- 17
            data_sub      <- data_set[ ( data_set$Time > time_divs[[j]][1] ) & ( data_set$Time < time_divs[[j]][2] ), ]
            data_sub_mu   <- mean( data_sub$Mean )
            # data_sub_sd   <- sqrt( sum( ( data_sub$Mean - data_sub$LB ) ** 2 ) / nrow( data_sub ) )
            data_sub_sd   <- sum( abs( data_sub$Mean - data_sub$LB ) * sqrt( green_data_N ) / 1.96 ) / nrow( data_sub )        # assumes Green's error bars are 95% CI, so +/- 1.96*SE
            # print( c( data_sub_mu, data_sub_sd ) )
            # print( seq( min( sims_sub_perc ), max( sims_sub_perc ), 0.01 ) )
            
            # ----- Stat #1: KS test statistic ----- #
            # data_sub_rand <- rtruncnorm( n=1000, a = 0, b = 100, mean = data_sub_mu, sd = data_sub_sd )
            
            # ks_res        <- ks.test( x = sims_sub_perc,
                                      # y = data_sub_rand,
                                      # alternative = "t", exact = FALSE )
            
            # # ----- Stat #2: KS test statistic normalized ----- #
            # ks_stat_norm  <- ks_res$statistic / ( data_sub_sd )
            
            # ----- Stat #3: Z-score ----- #
            # alt_stat      <- abs( sims_sub_perc_mu - data_sub_mu ) / ( data_sub_sd / sqrt( green_data_N ) )         # Z-score
            alt_stat      <- abs( sims_sub_perc_mu - data_sub_mu ) / ( data_sub_sd )         # Z-score
            # print( alt_stat )
            
            # ----- Stat #4: Z-score 2 ----- #
            
            z_score <- 0
            # print( data_sub )
            for( k in 1:nrow( data_sub ) ) {
                # print( data_sub$Time )
                pos_time_match <- which.min( abs( data_sub$Time[k] - sims_sub$time ) )
                z_score <- z_score + abs( data_sub$Mean[k] - sims_sub_perc[pos_time_match] ) / ( ( data_sub$Mean[k] - data_sub$LB[k] ) * sqrt( green_data_N ) )        # define a score as the average absolute value of z-score
                # z_score <- z_score / nrow( data_sub )
            }
            z_score <- z_score / nrow( data_sub )
            
        }
        
        # ----- Stat #1: KS test statistic ----- #
        # stat_sum      <- stat_sum + ks_res$statistic
        
        # # ----- Stat #2: KS test statistic normalized ----- #
        # stat_sum_norm <- stat_sum_norm + ks_stat_norm
        
        # ----- Stat #3: Z-score ----- #
        stat_sum_alt  <- stat_sum_alt + alt_stat
        
        # ----- Stat #4: Z-score 2 ----- #
        stat_sum_z    <- stat_sum_z + z_score
        
        # if( tolower( data_var ) == "cfu" ) {
            # cfu_stat_sum      <- stat_sum + ks_res$statistic
            # cfu_stat_sum_norm <- cfu_stat_sum_norm + ks_stat_norm
            # cfu_stat_sum_alt  <- cfu_stat_sum_alt + alt_stat
        # } else if( tolower( data_var ) == "tregpop" ) {
            # treg_stat_sum      <- stat_sum + ks_res$statistic
            # treg_stat_sum_norm <- treg_stat_sum_norm + ks_stat_norm
            # treg_stat_sum_alt  <- treg_stat_sum_alt + alt_stat
        # }
        
        
    }
    
    # return( stat_sum )
    # return( c( stat_sum, stat_sum_norm, stat_sum_alt, stat_sum_z ) )
    # return( c( stat_sum, stat_sum_alt, stat_sum_z ) )
    return( c( stat_sum_alt, stat_sum_z ) )
    # return( c( stat_sum, stat_sum_norm, stat_sum_alt, cfu_stat_sum, cfu_stat_sum_norm, cfu_stat_sum_alt, treg_stat_sum, treg_stat_sum_norm, treg_stat_sum_alt ) )
    
}


# flynn_stat_list <- flynn_stat_list_norm <- flynn_stat_list_alt <- flynn_stat_list_z <- list()
# green_stat_list <- green_stat_list_norm <- green_stat_list_alt <- green_stat_list_z <- list()
flynn_stat_list <- flynn_stat_list_alt <- flynn_stat_list_z <- list()
green_stat_list <- green_stat_list_alt <- green_stat_list_z <- list()

# for ( i in seq_along( sims[1:10] ) ) {
for ( i in seq_along( sims ) ) {
    
    # print( i )
    if ( i %% 100 == 0 ) {
        print( paste( i, " of ", length(sims), sep="" ) )
    }
    
    sims_df     <- data.frame( sims[[i]] )
    # print( sims_df )
    
    if ( i %in% dead_sims ) {
        
        # flynn_stat_list     <- c( flynn_stat_list, NA )
        # green_stat_list     <- c( green_stat_list, NA )
        
        flynn_stat_list_alt <- c( flynn_stat_list_alt, NA )
        green_stat_list_alt <- c( green_stat_list_alt, NA )
        
        flynn_stat_list_z <- c( flynn_stat_list_z, NA )
        green_stat_list_z <- c( green_stat_list_z, NA )
        
    } else {
        
        flynn_stat_sum  <- calc_stat_sum( sims_df, flynn_time_divs, flynn_data, "cfu" )
        green_stat_sum  <- calc_stat_sum( sims_df, green_time_divs, green_data, "tregpop" )
        # print( flynn_stat_sum )
        
        # flynn_stat_list     <- c( flynn_stat_list, flynn_stat_sum[1] )
        # green_stat_list     <- c( green_stat_list, green_stat_sum[1] )
        
        # flynn_stat_list_norm <- c( flynn_stat_list_norm, flynn_stat_sum[2] )
        # green_stat_list_norm <- c( green_stat_list_norm, green_stat_sum[2] )
        
        flynn_stat_list_alt <- c( flynn_stat_list_alt, flynn_stat_sum[1] )
        green_stat_list_alt <- c( green_stat_list_alt, green_stat_sum[1] )
        
        flynn_stat_list_z <- c( flynn_stat_list_z, flynn_stat_sum[2] )
        green_stat_list_z <- c( green_stat_list_z, green_stat_sum[2] )
    
    }
}


# cutoff      <- quantile( as.numeric( flynn_stat_list ), 0.01 )
# accept_pos  <- as.vector( which( flynn_stat_list < cutoff ) )
# cutoff      <- quantile( as.numeric( green_stat_list ), 0.01 )
# accept_pos  <- as.vector( which( green_stat_list < cutoff ) )


# ----- Stat #1: KS test statistic ----- #

# combo_score <- as.numeric( green_stat_list ) + as.numeric( flynn_stat_list )
# cutoff      <- quantile( as.numeric( combo_score ), accept_rate, na.rm=TRUE )
# accept_pos  <- as.vector( which( combo_score < cutoff ) )

# posts       <- priors[ accept_pos, ]
# posts       <- data.frame( posts, position=accept_pos, combo_score=combo_score[ accept_pos ] )

# sims_accept <- lapply( accept_pos, function(x){ data.frame( run_num=x, data.frame( sims[x] )[ c( "time", "mtb_act_gra_t", "imm_nai_t", "imm_pro_t", "imm_ant_t" ) ] ) } )
# sims_accept <- rbindlist( sims_accept )
# sims_accept$log_mtb_act_gra_t <- log10( sims_accept$mtb_act_gra_t + 1 )
# sims_accept$t_cell_perc       <- 100 * sims_accept$imm_ant_t / ( sims_accept$imm_nai_t + sims_accept$imm_pro_t + sims_accept$imm_ant_t )


# # ----- Stat #2: KS test statistic normalized ----- #

# combo_score_norm <- as.numeric( green_stat_list_norm ) + as.numeric( flynn_stat_list_norm )
# cutoff_norm      <- quantile( as.numeric( combo_score_norm ), 0.01, na.rm=TRUE )
# accept_pos_norm  <- as.vector( which( combo_score_norm < cutoff_norm ) )

# posts_norm       <- priors[ accept_pos_norm, ]
# posts_norm       <- data.frame( posts_norm, position=accept_pos_norm, combo_score=combo_score_norm[ accept_pos_norm ] )

# sims_accept_norm <- lapply( accept_pos_norm, function(x){ data.frame( run_num=x, data.frame( sims[x] )[ c( "time", "mtb_act_gra_t", "imm_nai_t", "imm_pro_t", "imm_ant_t" ) ] ) } )
# sims_accept_norm <- rbindlist( sims_accept_norm )
# sims_accept_norm$log_mtb_act_gra_t <- log10( sims_accept_norm$mtb_act_gra_t + 1 )
# sims_accept_norm$t_cell_perc       <- 100 * sims_accept_norm$imm_ant_t / ( sims_accept_norm$imm_nai_t + sims_accept_norm$imm_pro_t + sims_accept_norm$imm_ant_t )


# ----- Stat #3: alt score ----- #

combo_score_alt <- as.numeric( green_stat_list_alt ) + as.numeric( flynn_stat_list_alt )
cutoff_alt      <- quantile( as.numeric( combo_score_alt ), accept_rate, na.rm=TRUE )
accept_pos_alt  <- as.vector( which( combo_score_alt < cutoff_alt ) )

posts_alt   <- priors[ accept_pos_alt, ]
posts_alt   <- data.frame( posts_alt, position=accept_pos_alt, combo_score=combo_score_alt[ accept_pos_alt ] )

sims_accept_alt <- lapply( accept_pos_alt, function(x){ data.frame( run_num=x, data.frame( sims[x] )[ c( "time", "mtb_act_gra_t", "imm_nai_t", "imm_pro_t", "imm_ant_t" ) ] ) } )
sims_accept_alt <- rbindlist( sims_accept_alt )
sims_accept_alt$log_mtb_act_gra_t <- log10( sims_accept_alt$mtb_act_gra_t + 1 )
sims_accept_alt$t_cell_perc       <- 100 * sims_accept_alt$imm_ant_t / ( sims_accept_alt$imm_nai_t + sims_accept_alt$imm_pro_t + sims_accept_alt$imm_ant_t )


# ----- Stat #4: Z-score ----- #

combo_score_z <- 1 * as.numeric( green_stat_list_z ) + 1 * as.numeric( flynn_stat_list_z )
cutoff_z      <- quantile( as.numeric( combo_score_z ), accept_rate, na.rm=TRUE )
accept_pos_z  <- as.vector( which( combo_score_z <= cutoff_z ) )

posts_z   <- priors[ accept_pos_z, ]
posts_z   <- data.frame( posts_z, position=accept_pos_z, combo_score=combo_score_z[ accept_pos_z ] )

sims_accept_z <- lapply( accept_pos_z, function(x){ data.frame( run_num=x, data.frame( sims[x] )[ c( "time", "mtb_act_gra_t", "imm_nai_t", "imm_pro_t", "imm_ant_t" ) ] ) } )
sims_accept_z <- rbindlist( sims_accept_z )
sims_accept_z$log_mtb_act_gra_t <- log10( sims_accept_z$mtb_act_gra_t + 1 )
sims_accept_z$t_cell_perc       <- 100 * sims_accept_z$imm_ant_t / ( sims_accept_z$imm_nai_t + sims_accept_z$imm_pro_t + sims_accept_z$imm_ant_t )


# ----- Additional stats: Flynn data fit alone ----- #

# flynn_score   <- as.numeric( flynn_stat_list )
# flynn_min_pos <- as.vector( which.min( flynn_score ) )

# flynn_score_norm   <- as.numeric( flynn_stat_list_norm )
# flynn_min_pos_norm <- as.vector( which.min( flynn_score_norm ) )

flynn_score_alt   <- as.numeric( flynn_stat_list_alt )
flynn_min_pos_alt <- as.vector( which.min( flynn_score_alt ) )

flynn_score_z   <- as.numeric( flynn_stat_list_z )
flynn_min_pos_z <- as.vector( which.min( flynn_score_z ) )


# ----- Additional stats: Green data fit alone ----- #

# green_score   <- as.numeric( green_stat_list )
# green_min_pos <- as.vector( which.min( green_score ) )

# green_score_norm   <- as.numeric( green_stat_list_norm )
# green_min_pos_norm <- as.vector( which.min( green_score_norm ) )

green_score_alt   <- as.numeric( green_stat_list_alt )
green_min_pos_alt <- as.vector( which.min( green_score_alt ) )

green_score_z   <- as.numeric( green_stat_list_z )
green_min_pos_z <- as.vector( which.min( green_score_z ) )


# unique_time <- as.character( unique( sims_accept$time ) )
# log_mtb_act_gra_t_mean <- t_cell_perc_mean <- c()
# log_mtb_act_gra_t_mean_alt <- t_cell_perc_mean_alt <- c()
# for( i in 1:length( unique_time ) ) {
    
    # log_mtb_act_gra_t_mean     <- c( log_mtb_act_gra_t_mean,     mean( sims_accept[ which( as.character( sims_accept$time ) == unique_time[i] ), ]$log_mtb_act_gra_t ) )
    # log_mtb_act_gra_t_mean_alt <- c( log_mtb_act_gra_t_mean_alt, mean( sims_accept_alt[ which( as.character( sims_accept_alt$time ) == unique_time[i] ), ]$log_mtb_act_gra_t ) )
    
    # t_cell_perc_mean           <- c( t_cell_perc_mean,           mean( sims_accept[ which( as.character( sims_accept$time ) == unique_time[i] ), ]$t_cell_perc ) )
    # t_cell_perc_mean_alt       <- c( t_cell_perc_mean_alt,       mean( sims_accept_alt[ which( as.character( sims_accept_alt$time ) == unique_time[i] ), ]$t_cell_perc ) )
# }
# sims_accept_means     <- data.frame( time=as.numeric( as.character( unique_time ) ), log_mtb_act_gra_t_mean=log_mtb_act_gra_t_mean,     t_cell_perc_mean=t_cell_perc_mean )
# sims_accept_means_alt <- data.frame( time=as.numeric( as.character( unique_time ) ), log_mtb_act_gra_t_mean=log_mtb_act_gra_t_mean_alt, t_cell_perc_mean=t_cell_perc_mean_alt )


# p1 <- ggplot() +
    # theme_bw() +
    # theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
    # # geom_line( data=subset( sims_accept, time < round( max( flynn_data$Time ) / 1000 ) * 1000 ), aes( x=time, y=log_mtb_act_gra_t, group=run_num, colour=factor(run_num) ) ) +
    # geom_vline( xintercept=unique( unlist( flynn_time_divs ) ), linetype="dashed", size=0.75, colour="gray" ) +
    # # geom_line( data=subset( sims_accept, time < 750 ), aes( x=time, y=log_mtb_act_gra_t, group=run_num, colour=factor(run_num) ) ) +
    # geom_line( data=subset( sims_accept, time < 750 ), aes( x=time, y=log_mtb_act_gra_t, group=run_num ), colour="gray60", alpha=0.5 ) +
    # geom_point( data=data.frame( Time=flynn_data$Time, log_CFU=log10( flynn_data$CFU ) ), aes( x=Time, y=log_CFU ) ) +
    # # geom_line( data=subset( sims_accept_means, time < 750 ), aes( x=time, y=log_mtb_act_gra_t_mean ), colour="blue", size=0.75 ) +
    # # geom_line( data=subset( sims_accept, ( run_num == which.min( as.numeric( combo_score ) ) ) & ( time < 750 ) ), aes( x=time, y=log_mtb_act_gra_t ), colour="red", size=0.75 ) +
    # geom_line( data=subset( sims_accept, ( run_num == accept_pos[ which.min( combo_score[ accept_pos ] ) ] ) & ( time < 750 ) ), aes( x=time, y=log_mtb_act_gra_t ), colour="blue", size=1.5, alpha=0.9 ) +
    # # geom_line( data=subset( sims_accept, ( run_num == accept_pos[ which.max( combo_score[ accept_pos ] ) ] ) & ( time < 750 ) ), aes( x=time, y=log_mtb_act_gra_t ), colour="black", size=0.75 ) +
    # # geom_line( data=subset( sims_accept, ( run_num == flynn_min_pos ) & ( time < 750 ) ), aes( x=time, y=log_mtb_act_gra_t ), colour="blue", size=0.75 ) +
    # xlab( "Time (days)" ) +
    # ylab( "log10 CFU" ) +
    # ylim( -0.1, log10( mtb_max_gra ) )

# p2 <- ggplot() +
    # theme_bw() +
    # theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
    # # geom_line( data=subset( sims_accept, time < round( max( green_data$Time ) / 100 ) * 100 ), aes( x=time, y=t_cell_perc, group=run_num, colour=factor(run_num) ) ) +
    # geom_vline( xintercept=unique( unlist( green_time_divs ) ), linetype="dashed", size=0.75, colour="gray" ) +
    # # geom_line( data=subset( sims_accept, time < 200 ), aes( x=time, y=t_cell_perc, group=run_num, colour=factor(run_num) ) ) +
    # geom_line( data=subset( sims_accept, time < 200 ), aes( x=time, y=t_cell_perc, group=run_num ), colour="gray60", alpha=0.5 ) +
    # geom_point( data=data.frame( Time=green_data$Time, Treg_perc=green_data$Mean ), aes( x=Time, y=Treg_perc ) ) +
    # # geom_line( data=subset( sims_accept_means, time < 200 ), aes( x=time, y=t_cell_perc_mean ), colour="blue", size=0.75 ) +
    # # geom_line( data=subset( sims_accept, ( run_num == which.min( as.numeric( combo_score ) ) ) & ( time < 200 ) ), aes( x=time, y=t_cell_perc ), colour="red", size=0.75 ) +
    # geom_line( data=subset( sims_accept, ( run_num == accept_pos[ which.min( combo_score[ accept_pos ] ) ] ) & ( time < 200 ) ), aes( x=time, y=t_cell_perc ), colour="blue", size=1.5, alpha=0.9 ) +
    # # geom_line( data=subset( sims_accept, ( run_num == accept_pos[ which.max( combo_score[ accept_pos ] ) ] ) & ( time < 200 ) ), aes( x=time, y=t_cell_perc ), colour="black", size=0.75 ) +
    # # geom_line( data=subset( sims_accept, ( run_num == green_min_pos ) & ( time < 200 ) ), aes( x=time, y=t_cell_perc ), colour="blue", size=0.75 ) +
    # geom_errorbar( data=data.frame( Time=green_data$Time, LB_SD=green_data$LB_SD, UB_SD=green_data$UB_SD ), aes( x=Time, ymin=LB_SD, ymax=UB_SD ), size=0.75, width=0, colour="black" ) +
    # # geom_errorbar( data=data.frame( Time=green_data$Time, LB=green_data$LB, UB=green_data$UB ), aes( x=Time, ymin=LB, ymax=UB ), size=0.75, width=4, colour="black" ) +
    # # ylim( -5, 100 ) +
    # xlab( "Time (days)" ) +
    # ylab( "Treg %age of CD4" )


p1_alt <- ggplot() +
    theme_bw() +
    theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
    # geom_line( data=subset( sims_accept, time < round( max( flynn_data$Time ) / 1000 ) * 1000 ), aes( x=time, y=log_mtb_act_gra_t, group=run_num, colour=factor(run_num) ) ) +
    geom_vline( xintercept=unique( unlist( flynn_time_divs ) ), linetype="dashed", size=0.75, colour="gray" ) +
    # geom_line( data=subset( sims_accept_alt, time < 750 ), aes( x=time, y=log_mtb_act_gra_t, group=run_num, colour=factor(run_num) ) ) +
    geom_line( data=subset( sims_accept_alt, time < 750 ), aes( x=time, y=log_mtb_act_gra_t, group=run_num ), colour="gray60", alpha=0.5 ) +
    geom_point( data=data.frame( Time=flynn_data$Time, log_CFU=log10( flynn_data$CFU ) ), aes( x=Time, y=log_CFU ) ) +
    # geom_line( data=subset( sims_accept_means_alt, time < 750 ), aes( x=time, y=log_mtb_act_gra_t_mean ), colour="blue", size=0.75 ) +
    # geom_line( data=subset( sims_accept_alt, ( run_num == which.min( as.numeric( combo_score_alt ) ) ) & ( time < 750 ) ), aes( x=time, y=log_mtb_act_gra_t ), colour="red", size=0.75 ) +
    geom_line( data=subset( sims_accept_alt, ( run_num == accept_pos_alt[ which.min( combo_score_alt[ accept_pos_alt ] ) ] ) & ( time < 750 ) ), aes( x=time, y=log_mtb_act_gra_t ), colour="blue", size=1.5, alpha=0.9 ) +
    # geom_line( data=subset( sims_accept_alt, ( run_num == accept_pos_alt[ which.max( combo_score_alt[ accept_pos_alt ] ) ] ) & ( time < 750 ) ), aes( x=time, y=log_mtb_act_gra_t ), colour="black", size=0.75 ) +
    # geom_line( data=subset( sims_accept_alt, ( run_num == flynn_min_pos_alt ) & ( time < 750 ) ), aes( x=time, y=log_mtb_act_gra_t ), colour="blue", size=0.75 ) +
    xlab( "Time (days)" ) +
    ylab( "log10 CFU" ) +
    ylim( -0.1, log10( mtb_max_gra ) )

p2_alt <- ggplot() +
    theme_bw() +
    theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
    # geom_line( data=subset( sims_accept, time < round( max( green_data$Time ) / 100 ) * 100 ), aes( x=time, y=t_cell_perc, group=run_num, colour=factor(run_num) ) ) +
    geom_vline( xintercept=unique( unlist( green_time_divs ) ), linetype="dashed", size=0.75, colour="gray" ) +
    # geom_line( data=subset( sims_accept_alt, time < 200 ), aes( x=time, y=t_cell_perc, group=run_num, colour=factor(run_num) ) ) +
    geom_line( data=subset( sims_accept_alt, time < 200 ), aes( x=time, y=t_cell_perc, group=run_num ), colour="gray60", alpha=0.5 ) +
    geom_point( data=data.frame( Time=green_data$Time, Treg_perc=green_data$Mean ), aes( x=Time, y=Treg_perc ) ) +
    # geom_line( data=subset( sims_accept_means_alt, time < 200 ), aes( x=time, y=t_cell_perc_mean ), colour="blue", size=0.75 ) +
    # geom_line( data=subset( sims_accept_alt, ( run_num == which.min( as.numeric( combo_score_alt ) ) ) & ( time < 200 ) ), aes( x=time, y=t_cell_perc ), colour="red", size=0.75 ) +
    geom_line( data=subset( sims_accept_alt, ( run_num == accept_pos_alt[ which.min( combo_score_alt[ accept_pos_alt ] ) ] ) & ( time < 200 ) ), aes( x=time, y=t_cell_perc ), colour="blue", size=1.5, alpha=0.9 ) +
    # geom_line( data=subset( sims_accept_alt, ( run_num == accept_pos_alt[ which.max( combo_score_alt[ accept_pos_alt ] ) ] ) & ( time < 200 ) ), aes( x=time, y=t_cell_perc ), colour="black", size=0.75 ) +
    # geom_line( data=subset( sims_accept_alt, ( run_num == green_min_pos_alt ) & ( time < 200 ) ), aes( x=time, y=t_cell_perc ), colour="blue", size=0.75 ) +
    geom_errorbar( data=data.frame( Time=green_data$Time, LB_SD=green_data$LB_SD, UB_SD=green_data$UB_SD ), aes( x=Time, ymin=LB_SD, ymax=UB_SD ), size=0.75, width=0, colour="black" ) +
    # geom_errorbar( data=data.frame( Time=green_data$Time, LB=green_data$LB, UB=green_data$UB ), aes( x=Time, ymin=LB, ymax=UB ), size=0.75, width=4, colour="black" ) +
    # ylim( -5, 100 ) +
    xlab( "Time (days)" ) +
    ylab( "Treg %age of CD4" )


p1_z <- ggplot() +
    theme_bw() +
    theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
    # geom_line( data=subset( sims_accept, time < round( max( flynn_data$Time ) / 1000 ) * 1000 ), aes( x=time, y=log_mtb_act_gra_t, group=run_num, colour=factor(run_num) ) ) +
    geom_vline( xintercept=unique( unlist( flynn_time_divs ) ), linetype="dashed", size=0.75, colour="gray" ) +
    # geom_line( data=subset( sims_accept_z, time < 750 ), aes( x=time, y=log_mtb_act_gra_t, group=run_num, colour=factor(run_num) ) ) +
    geom_line( data=subset( sims_accept_z, time < 750 ), aes( x=time, y=log_mtb_act_gra_t, group=run_num ), colour="gray60", alpha=0.5 ) +
    geom_point( data=data.frame( Time=flynn_data$Time, log_CFU=log10( flynn_data$CFU ) ), aes( x=Time, y=log_CFU ) ) +
    # geom_line( data=subset( sims_accept_means_z, time < 750 ), aes( x=time, y=log_mtb_act_gra_t_mean ), colour="blue", size=0.75 ) +
    # geom_line( data=subset( sims_accept_z, ( run_num == which.min( as.numeric( combo_score_z ) ) ) & ( time < 750 ) ), aes( x=time, y=log_mtb_act_gra_t ), colour="red", size=0.75 ) +
    geom_line( data=subset( sims_accept_z, ( run_num == accept_pos_z[ which.min( combo_score_z[ accept_pos_z ] ) ] ) & ( time < 750 ) ), aes( x=time, y=log_mtb_act_gra_t ), colour="blue", size=1.5, alpha=0.9 ) +
    # geom_line( data=subset( sims_accept_z, ( run_num == accept_pos_z[ which.max( combo_score_z[ accept_pos_z ] ) ] ) & ( time < 750 ) ), aes( x=time, y=log_mtb_act_gra_t ), colour="black", size=0.75 ) +
    xlab( "Time (days)" ) +
    ylab( "log10 CFU" ) +
    ylim( -0.1, log10( mtb_max_gra ) )

p2_z <- ggplot() +
    theme_bw() +
    theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
    # geom_line( data=subset( sims_accept, time < round( max( green_data$Time ) / 100 ) * 100 ), aes( x=time, y=t_cell_perc, group=run_num, colour=factor(run_num) ) ) +
    geom_vline( xintercept=unique( unlist( green_time_divs ) ), linetype="dashed", size=0.75, colour="gray" ) +
    # geom_line( data=subset( sims_accept_z, time < 200 ), aes( x=time, y=t_cell_perc, group=run_num, colour=factor(run_num) ) ) +
    geom_line( data=subset( sims_accept_z, time < 200 ), aes( x=time, y=t_cell_perc, group=run_num ), colour="gray60", alpha=0.5 ) +
    geom_point( data=data.frame( Time=green_data$Time, Treg_perc=green_data$Mean ), aes( x=Time, y=Treg_perc ) ) +
    # geom_errorbar( data=data.frame( Time=green_data$Time, LB_SD=green_data$LB_SD, UB_SD=green_data$UB_SD ), aes( x=Time, ymin=LB_SD, ymax=UB_SD ), size=0.75, width=4, colour="gray60" ) +
    # geom_errorbar( data=data.frame( Time=green_data$Time, LB=green_data$LB, UB=green_data$UB ), aes( x=Time, ymin=LB, ymax=UB ), size=0.75, width=4, colour="black" ) +
    # # geom_line( data=subset( sims_accept_means_z, time < 200 ), aes( x=time, y=t_cell_perc_mean ), colour="blue", size=0.75 ) +
    # # geom_line( data=subset( sims_accept_z, ( run_num == which.min( as.numeric( combo_score_z ) ) ) & ( time < 200 ) ), aes( x=time, y=t_cell_perc ), colour="red", size=0.75 ) +
    geom_line( data=subset( sims_accept_z, ( run_num == accept_pos_z[ which.min( combo_score_z[ accept_pos_z ] ) ] ) & ( time < 200 ) ), aes( x=time, y=t_cell_perc ), colour="blue", size=1.5, alpha=0.9 ) +
    geom_errorbar( data=data.frame( Time=green_data$Time, LB_SD=green_data$LB_SD, UB_SD=green_data$UB_SD ), aes( x=Time, ymin=LB_SD, ymax=UB_SD ), size=0.75, width=0, colour="black" ) +
    # geom_line( data=subset( sims_accept_z, ( run_num == accept_pos_z[ which.max( combo_score_z[ accept_pos_z ] ) ] ) & ( time < 200 ) ), aes( x=time, y=t_cell_perc ), colour="black", size=0.75 ) +
    xlab( "Time (days)" ) +
    ylab( "Treg %age of CD4" )


# hist1d_1 <- ggplot() +
            # geom_histogram( data=data.frame( posts ), aes( x=cmi_kill_rate_log, y=..density.. ), color="black", fill="white", binwidth=1 ) +
            # geom_density( data=data.frame( posts ), aes( x=cmi_kill_rate_log, y=..density.. ), color="blue" ) +
            # xlim( round( min( data.frame( posts )$cmi_kill_rate_log ) ), round( max( data.frame( posts )$cmi_kill_rate_log ) ) ) +
            # ggtitle( "Posterior distribution" ) +
            # theme_bw()
# hist1d_2 <- ggplot() +
            # geom_histogram( data=data.frame( posts ), aes( x=nai_decay_rate_log, y=..density.. ), color="black", fill="white", binwidth=1 ) +
            # geom_density( data=data.frame( posts ), aes( x=nai_decay_rate_log, y=..density.. ), color="blue" ) +
            # xlim( round( min( data.frame( posts )$nai_decay_rate_log ) ), round( max( data.frame( posts )$nai_decay_rate_log ) ) ) +
            # ggtitle( "Posterior distribution" ) +
            # theme_bw()


# ################################################
# # ----- Plot #1a: KS test statistic fits ----- #
# ################################################

# pdf( height=5, width=8, "output_fitting.pdf" )
# print( p1 )
# print( p2 )
# dev.off()


# # ----- Plot #1b: KS test statistic param dists ----- #

# hist2d_1  <- kde2d( x=posts[,1], y=posts[,3] )
# hist2d_1b <- kde2d( x=posts[,2], y=posts[,3] )
# hist2d_2  <- kde2d( x=posts[,4], y=posts[,5] )
# hist2d_2b <- kde2d( x=posts[,3], y=posts[,6] )
# # hist2d_3  <- kde2d( x=posts[,7], y=posts[,8] )

init_rates_df <- data.frame( param = names(init_rates), init_rates )

# # init_rates_df$param <- as.character( init_rates_df$param )
# # init_rates_df$param[ init_rates_df$param=="growth_rate_log" ]       <- "growth_rate_log"
# # init_rates_df$param[ init_rates_df$param=="cmi_kill_rate_log" ]     <- "killing_rate_log"
# # init_rates_df$param[ init_rates_df$param=="imm_source_rate_log" ]   <- "imm_priming_rate_log"
# # init_rates_df$param[ init_rates_df$param=="nai2pro_diff_rate_log" ] <- "Teff_diff_rate_log"
# # init_rates_df$param[ init_rates_df$param=="nai2ant_diff_rate_log" ] <- "Treg_diff_rate_log"
# # init_rates_df$param[ init_rates_df$param=="nai_decay_rate_log" ]    <- "Tpri_decay_rate_log"
# # init_rates_df$param[ init_rates_df$param=="pro_decay_rate_log" ]    <- "Teff_decay_rate_log"
# # init_rates_df$param[ init_rates_df$param=="ant_decay_rate_log" ]    <- "Treg_decay_rate_log"


# posts_melt <- melt( posts[,1:(num_params+1)], id="position" )
# colnames( posts_melt ) <- c( "position", "param", "value" )

# # posts_melt$param <- as.character( posts_melt$param )
# # posts_melt$param[ posts_melt$param=="growth_rate_log" ]       <- "growth_rate_log"
# # posts_melt$param[ posts_melt$param=="cmi_kill_rate_log" ]     <- "killing_rate_log"
# # posts_melt$param[ posts_melt$param=="imm_source_rate_log" ]   <- "imm_priming_rate_log"
# # posts_melt$param[ posts_melt$param=="nai2pro_diff_rate_log" ] <- "Teff_diff_rate_log"
# # posts_melt$param[ posts_melt$param=="nai2ant_diff_rate_log" ] <- "Treg_diff_rate_log"
# # posts_melt$param[ posts_melt$param=="nai_decay_rate_log" ]    <- "Tpri_decay_rate_log"
# # posts_melt$param[ posts_melt$param=="pro_decay_rate_log" ]    <- "Teff_decay_rate_log"
# # posts_melt$param[ posts_melt$param=="ant_decay_rate_log" ]    <- "Treg_decay_rate_log"

# cdat <- ddply( posts_melt, "param", summarise, mean=mean( value ), median=median( value ) ) # , bin=diff( range( value ) ) )


# hist1d_all <- ggplot() +
    # theme_bw() +
    # theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
    # geom_histogram( data=posts_melt, aes( x=value, y=..density.. ), colour="gray", fill="gray" ) + # , binwidth=diff(range(c(1,2)))/10 ) +
    # geom_density( data=posts_melt, aes( x=value, y=..density.. ) ) +
    # geom_vline( data=cdat, aes( xintercept=mean ), linetype="solid", size=2, alpha=0.5 ) +
    # # geom_vline( data=cdat, aes( xintercept=median ), linetype="dotted", size=1 ) +
    # geom_vline( data=init_rates_df, aes( xintercept=init_rates ), colour="green", size=2, alpha=0.7 ) +
    # geom_vline( data=subset( posts_melt, position==accept_pos[ which.min( combo_score[ accept_pos ] ) ] ), aes( xintercept=value ), colour="blue", size=2, alpha=0.5 ) +
    # facet_wrap( ~param, scales="free" )

# hist_list <- list()
# col_names <- colnames( posts )[1:num_params]
# highest <- subset( posts_melt, position==accept_pos[ which.min( combo_score[ accept_pos ] ) ] ) 

# pdf( height=5, width=8, "output_fitting_hists.pdf" )

# print( hist1d_all )

# for( i in 1:ncol( posts[,1:num_params] ) ) {
    # # print( i )
    # # foo <- qplot( posts[,i], geom="histogram" )
    # foo <- ggplot() +
        # theme_bw() +
        # theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
        # geom_histogram( data=data.frame( value=posts[,i] ), aes( x=value, y=..density.. ), colour="gray", fill="gray", binwidth=diff( range( posts[,i] ) )/10 ) +       # histogram
        # geom_density( data=data.frame( value=posts[,i] ), aes( x=value, y=..density.. ) ) +     # kernel density
        # geom_vline( data=cdat[i,], aes( xintercept=mean ), linetype="solid", size=2, alpha=0.5 ) +        # mean of histogram
        # geom_vline( data=init_rates_df[i,], aes( xintercept=init_rates ), colour="green", size=2, alpha=0.7 ) +
        # geom_vline( data=highest[i,], aes( xintercept=value ), colour="blue", size=2, alpha=0.5 ) +         # value for best overall fitting
        # xlab( colnames( posts )[i] )
    # print( foo )
# }

# contour( hist2d_1, xlab=dimnames(posts)[[2]][1], ylab=dimnames(posts)[[2]][3] )
# # # print( hist1d_1 )
# points( x=cdat$mean[1], y=cdat$mean[3], col="black", pch=16 )
# points( x=init_rates_df[1,2], y=init_rates_df[3,2], col="green", pch=16 )
# points( x=highest$value[1], y=highest$value[3], col="blue", pch=16 )


# contour( hist2d_1b, xlab=dimnames(posts)[[2]][2], ylab=dimnames(posts)[[2]][3] )
# points( x=cdat$mean[2], y=cdat$mean[3], col="black", pch=16 )
# points( x=init_rates_df[2,2], y=init_rates_df[3,2], col="green", pch=16 )
# points( x=highest$value[2], y=highest$value[3], col="blue", pch=16 )

# # hist( data.frame(posts)[ , grep( "cmi_kill_rate_log", colnames( posts ) ) ], prob=T, main="", xlab="cmi_kill_rate_log" ); lines( density( posts[,2] ) )
# contour( hist2d_2, xlab=dimnames(posts)[[2]][4], ylab=dimnames(posts)[[2]][5] )
# # # print( hist1d_2 )
# points( x=cdat$mean[4], y=cdat$mean[5], col="black", pch=16 )
# points( x=init_rates_df[4,2], y=init_rates_df[5,2], col="green", pch=16 )
# points( x=highest$value[4], y=highest$value[5], col="blue", pch=16 )

# contour( hist2d_2b, xlab=dimnames(posts)[[2]][3], ylab=dimnames(posts)[[2]][6] )
# points( x=cdat$mean[3], y=cdat$mean[6], col="black", pch=16 )
# points( x=init_rates_df[3,2], y=init_rates_df[6,2], col="green", pch=16 )
# points( x=highest$value[3], y=highest$value[6], col="blue", pch=16 )

# # # hist( data.frame(posts)[ , grep( "nai_decay_rate_log", colnames( posts ) ) ], prob=T, main="", xlab="nai_decay_rate_log" ); lines( density( posts[,6] ) )
# # contour( hist2d_3, xlab=dimnames(posts)[[2]][7], ylab=dimnames(posts)[[2]][8] )
# # points( x=cdat$mean[7], y=cdat$mean[8], col="black", pch=16 )
# # points( x=init_rates_df[7,2], y=init_rates_df[8,2], col="green", pch=16 )
# # points( x=highest$value[7], y=highest$value[8], col="blue", pch=16 )
# dev.off()


# # ----- Output #1c: KS test statistic posterior distributions ----- #

# # write.csv( data.frame( position=accept_pos, posts, combo_score=combo_score[ accept_pos ] ), file="posts.csv" )
# write.csv( posts, file="posts.csv" )
# write.csv( subset( sims_accept, time < 750 ), file="time_courses.csv", row.names=F )

# # write.csv( posts_z, file="posts_z.csv" )


########################################
# ----- Plot #2a: alt score fits ----- #
########################################

pdf( height=5, width=8, "output_fitting_alt.pdf" )
print( p1_alt )
print( p2_alt )
dev.off()


# ----- Plot #2b: Z-score param dists ----- #

hist2d_1_alt  <- kde2d( x=posts_alt[,1], y=posts_alt[,3] )
hist2d_1b_alt <- kde2d( x=posts_alt[,2], y=posts_alt[,3] )
hist2d_2_alt  <- kde2d( x=posts_alt[,4], y=posts_alt[,5] )
hist2d_2b_alt <- kde2d( x=posts_alt[,3], y=posts_alt[,6] )
# hist2d_3_alt  <- kde2d( x=posts_alt[,7], y=posts_alt[,8] )

posts_alt_melt <- melt( posts_alt[,1:(num_params+1)], id="position" )
colnames( posts_alt_melt ) <- c( "position", "param", "value" )
cdat <- ddply( posts_alt_melt, "param", summarise, mean=mean( value ), median=median( value ) ) # , bin=diff( range( value ) ) )

hist1d_all_alt <- ggplot() +
    theme_bw() +
    theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
    geom_histogram( data=posts_alt_melt, aes( x=value, y=..density.. ), colour="gray", fill="gray" ) + # , binwidth=diff(range(c(1,2)))/10 ) +
    geom_density( data=posts_alt_melt, aes( x=value, y=..density.. ) ) +
    geom_vline( data=cdat, aes( xintercept=mean ), linetype="solid", size=2, alpha=0.5 ) +
    # geom_vline( data=cdat, aes( xintercept=median ), linetype="dotted", size=1 ) +
    geom_vline( data=init_rates_df, aes( xintercept=init_rates ), colour="green", size=2, alpha=0.7 ) +
    geom_vline( data=subset( posts_alt_melt, position==accept_pos_alt[ which.min( combo_score_alt[ accept_pos_alt ] ) ] ), aes( xintercept=value ), colour="blue", size=2, alpha=0.5 ) +
    facet_wrap( ~param, scales="free" )

hist_list <- list()
col_names <- colnames( posts_alt )[1:num_params]
highest <- subset( posts_alt_melt, position==accept_pos_alt[ which.min( combo_score_alt[ accept_pos_alt ] ) ] )


pdf( height=5, width=8, "output_fitting_alt_hists.pdf" )

print( hist1d_all_alt )

for( i in 1:ncol( posts_alt[,1:num_params] ) ) {
    # print( i )
    # foo <- qplot( posts_alt[,i], geom="histogram" )
    foo <- ggplot() +
        theme_bw() +
        theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
        geom_histogram( data=data.frame( value=posts_alt[,i] ), aes( x=value, y=..density.. ), colour="gray", fill="gray", binwidth=diff( range( posts_alt[,i] ) )/10 ) +
        geom_density( data=data.frame( value=posts_alt[,i] ), aes( x=value, y=..density.. ) ) +
        geom_vline( data=cdat[i,], aes( xintercept=mean ), linetype="solid", size=2, alpha=0.5 ) +
        geom_vline( data=init_rates_df[i,], aes( xintercept=init_rates ), colour="green", size=2, alpha=0.7 ) +
        geom_vline( data=highest[i,], aes( xintercept=value ), colour="blue", size=2, alpha=0.5 ) +
        xlab( colnames( posts_alt )[i] )
    print( foo )
}

contour( hist2d_1_alt, xlab=dimnames(posts_alt)[[2]][1], ylab=dimnames(posts_alt)[[2]][3] )
# # print( hist1d_1 )
points( x=cdat$mean[1], y=cdat$mean[3], col="black", pch=16 )
points( x=init_rates_df[1,2], y=init_rates_df[3,2], col="green", pch=16 )
points( x=highest$value[1], y=highest$value[3], col="blue", pch=16 )

contour( hist2d_1b_alt, xlab=dimnames(posts_alt)[[2]][2], ylab=dimnames(posts_alt)[[2]][3] )
points( x=cdat$mean[2], y=cdat$mean[3], col="black", pch=16 )
points( x=init_rates_df[2,2], y=init_rates_df[3,2], col="green", pch=16 )
points( x=highest$value[2], y=highest$value[3], col="blue", pch=16 )

# hist( data.frame(posts_alt)[ , grep( "cmi_kill_rate_log", colnames( posts_alt ) ) ], prob=T, main="", xlab="cmi_kill_rate_log" ); lines( density( posts_alt[,2] ) )
contour( hist2d_2_alt, xlab=dimnames(posts_alt)[[2]][4], ylab=dimnames(posts_alt)[[2]][5] )
# # print( hist1d_2 )
points( x=cdat$mean[4], y=cdat$mean[5], col="black", pch=16 )
points( x=init_rates_df[4,2], y=init_rates_df[5,2], col="green", pch=16 )
points( x=highest$value[4], y=highest$value[5], col="blue", pch=16 )

contour( hist2d_2b_alt, xlab=dimnames(posts_alt)[[2]][3], ylab=dimnames(posts_alt)[[2]][6] )
points( x=cdat$mean[3], y=cdat$mean[6], col="black", pch=16 )
points( x=init_rates_df[3,2], y=init_rates_df[6,2], col="green", pch=16 )
points( x=highest$value[3], y=highest$value[6], col="blue", pch=16 )

# # hist( data.frame(posts_alt)[ , grep( "nai_decay_rate_log", colnames( posts_alt ) ) ], prob=T, main="", xlab="nai_decay_rate_log" ); lines( density( posts_alt[,6] ) )
# contour( hist2d_3_alt, xlab=dimnames(posts_alt)[[2]][7], ylab=dimnames(posts_alt)[[2]][8] )
# points( x=cdat$mean[7], y=cdat$mean[8], col="black", pch=16 )
# points( x=init_rates_df[7,2], y=init_rates_df[8,2], col="green", pch=16 )
# points( x=highest$value[7], y=highest$value[8], col="blue", pch=16 )

dev.off()


# ----- Output #2c: Z-score posterior distributions ----- #

# write.csv( data.frame( position=accept_pos_alt, posts_alt, combo_score=combo_score_alt[ accept_pos_alt ] ), file="posts_alt.csv" )
write.csv( posts_alt, file="posts_alt.csv" )
write.csv( subset( sims_accept_alt, time < 750 ), file="time_courses_alt.csv", row.names=F )


#######################################
# ----- Plot #3a: Z-score2 fits ----- #
#######################################

pdf( height=5, width=8, "output_fitting_z.pdf" )
print( p1_z )
print( p2_z )
dev.off()


# ----- Plot #3b: Z-score2 param dists ----- #

hist2d_1_z  <- kde2d( x=posts_z[,1], y=posts_z[,3] )
hist2d_1b_z <- kde2d( x=posts_z[,2], y=posts_z[,3] )
hist2d_2_z  <- kde2d( x=posts_z[,4], y=posts_z[,5] )
hist2d_2b_z <- kde2d( x=posts_z[,3], y=posts_z[,6] )
# hist2d_3_z  <- kde2d( x=posts_z[,7], y=posts_z[,8] )

posts_z_melt <- melt( posts_z[,1:(num_params+1)], id="position" )
colnames( posts_z_melt ) <- c( "position", "param", "value" )
cdat <- ddply( posts_z_melt, "param", summarise, mean=mean( value ), median=median( value ) ) # , bin=diff( range( value ) ) )  # report cdat table as supplement

hist1d_all_z <- ggplot() +
    theme_bw() +
    theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
    geom_histogram( data=posts_z_melt, aes( x=value, y=..density.. ), colour="gray", fill="gray" ) + # , binwidth=diff(range(c(1,2)))/10 ) +
    geom_density( data=posts_z_melt, aes( x=value, y=..density.. ) ) +
    geom_vline( data=cdat, aes( xintercept=mean ), linetype="solid", colour="black", size=2, alpha=0.5 ) +
    # geom_vline( data=cdat, aes( xintercept=median ), linetype="dotted", size=1 ) +
    geom_vline( data=init_rates_df, aes( xintercept=init_rates ), colour="green", size=2, alpha=0.7 ) +
    geom_vline( data=subset( posts_z_melt, position==accept_pos_z[ which.min( combo_score_z[ accept_pos_z ] ) ] ), aes( xintercept=value ), colour="blue", size=2, alpha=0.5 ) +
    facet_wrap( ~param, scales="free" )

hist_list <- list()
col_names <- colnames( posts_z )[1:8]
highest <- subset( posts_z_melt, position==accept_pos_z[ which.min( combo_score_z[ accept_pos_z ] ) ] )
# report highest table as supplement


pdf( height=5, width=8, "output_fitting_z_hists.pdf" )

print( hist1d_all_z )

for( i in 1:ncol( posts_z[,1:num_params] ) ) {
    # print( i )
    # foo <- qplot( posts[,i], geom="histogram" )
    foo <- ggplot() +
        theme_bw() +
        theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
        geom_histogram( data=data.frame( value=posts_z[,i] ), aes( x=value, y=..density.. ), colour="gray", fill="gray", binwidth=diff( range( posts_z[,i] ) )/10 ) +       # histogram
        geom_density( data=data.frame( value=posts_z[,i] ), aes( x=value, y=..density.. ) ) +     # kernel density
        geom_vline( data=cdat[i,], aes( xintercept=mean ), linetype="solid", size=2, alpha=0.5 ) +        # mean of histogram
        geom_vline( data=init_rates_df[i,], aes( xintercept=init_rates ), colour="green", size=2, alpha=0.7 ) +
        geom_vline( data=highest[i,], aes( xintercept=value ), colour="blue", size=2, alpha=0.5 ) +         # value for best overall fitting
        xlab( colnames( posts_z )[i] )
    print( foo )
}

contour( hist2d_1_z, xlab=dimnames(posts_z)[[2]][1], ylab=dimnames(posts_z)[[2]][3] )
# # print( hist1d_1 )
points( x=cdat$mean[1], y=cdat$mean[3], col="black", pch=16 )
points( x=init_rates_df[1,2], y=init_rates_df[3,2], col="green", pch=16 )
points( x=highest$value[1], y=highest$value[3], col="blue", pch=16 )

contour( hist2d_1b_z, xlab=dimnames(posts_z)[[2]][2], ylab=dimnames(posts_z)[[2]][3] )
points( x=cdat$mean[2], y=cdat$mean[3], col="black", pch=16 )
points( x=init_rates_df[2,2], y=init_rates_df[3,2], col="green", pch=16 )
points( x=highest$value[2], y=highest$value[3], col="blue", pch=16 )

# hist( data.frame(posts)[ , grep( "cmi_kill_rate_log", colnames( posts ) ) ], prob=T, main="", xlab="cmi_kill_rate_log" ); lines( density( posts[,2] ) )
contour( hist2d_2_z, xlab=dimnames(posts_z)[[2]][4], ylab=dimnames(posts_z)[[2]][5] )
# # print( hist1d_2 )
points( x=cdat$mean[4], y=cdat$mean[5], col="black", pch=16 )
points( x=init_rates_df[4,2], y=init_rates_df[5,2], col="green", pch=16 )
points( x=highest$value[4], y=highest$value[5], col="blue", pch=16 )

contour( hist2d_2b_z, xlab=dimnames(posts_z)[[2]][3], ylab=dimnames(posts_z)[[2]][6] )
points( x=cdat$mean[3], y=cdat$mean[6], col="black", pch=16 )
points( x=init_rates_df[3,2], y=init_rates_df[6,2], col="green", pch=16 )
points( x=highest$value[3], y=highest$value[6], col="blue", pch=16 )

# # hist( data.frame(posts)[ , grep( "nai_decay_rate_log", colnames( posts ) ) ], prob=T, main="", xlab="nai_decay_rate_log" ); lines( density( posts[,6] ) )
# contour( hist2d_3_z, xlab=dimnames(posts)[[2]][7], ylab=dimnames(posts)[[2]][8] )
# points( x=cdat$mean[7], y=cdat$mean[8], col="black", pch=16 )
# points( x=init_rates_df[7,2], y=init_rates_df[8,2], col="green", pch=16 )
# points( x=highest$value[7], y=highest$value[8], col="blue", pch=16 )

dev.off()


# ----- Output #3c: Z-score2 posterior distributions ----- #

# write.csv( data.frame( position=accept_pos_alt, posts_alt, combo_score=combo_score_alt[ accept_pos_alt ] ), file="posts_alt.csv" )
write.csv( posts_z, file="posts_z.csv" )
write.csv( subset( sims_accept_z, time < 750 ), file="time_courses_z.csv", row.names=F )

#####################################################
#====================================================
# Downstream output: incidence
#====================================================
#####################################################

library(RColorBrewer)

# Take selected parameters and re-run simulations for longer time periods

start_time <- 0
end_time   <- 3650
intervals  <- 1001

timesteps  <- seq( from=start_time, to=end_time, length.out=intervals )

y0         <- c( mtb_qui_gra_t=mtb_qui_gra_init, mtb_act_gra_t=mtb_act_gra_init, imm_nai_t=imm_nai_init, imm_pro_t=imm_pro_init, imm_ant_t=imm_ant_init )

init_rates <- c(    growth_rate_log=growth_rate_log,
                    cmi_kill_rate_log=cmi_kill_rate_log,
                    imm_source_rate_log=imm_source_rate_log,
                    nai2pro_diff_rate_log=nai2pro_diff_rate_log,
                    nai2ant_diff_rate_log=nai2ant_diff_rate_log,
                    nai_decay_rate_log=nai_decay_rate_log,
                    pro_decay_rate_log=pro_decay_rate_log,
                    ant_decay_rate_log=ant_decay_rate_log )

# init_rates <- c( growth_rate_log=growth_rate_log, cmi_kill_rate_log=cmi_kill_rate_log, imm_source_rate_log=imm_source_rate_log, nai2pro_diff_rate_log=nai2pro_diff_rate_log, nai2ant_diff_rate_log=nai2ant_diff_rate_log, nai_decay_rate_log=nai_decay_rate_log, imm_ant_eff_log=imm_ant_eff_log )


# ####################################################
# ### KS #############################################
# ####################################################

# posts_list <- vector( mode="list", length=nrow( posts ) )
# for ( i in 1:length( posts_list ) ) {
    # # print( i )
    # posts_list[[ i ]] <- posts[ i, ]
# }

# # sims_posts <- llply( .data = posts_list, .fun = function( lhs_params ){
    # # init_rates[ ( names( lhs_params )[1:length(init_rates)] ) ] <- lhs_params[1:length(init_rates)]
    # # # init_rates[ names( lhs_params[ 1:length(init_rates) ] ) ] <- lhs_params
    # # init_rates <- unlist( init_rates )
    # # print( init_rates )
    # # ode( y=y0, times=timesteps, func=model_eq, parms=init_rates, rtol=1e-15, maxsteps=20000, method="lsoda" )
# # }, .progress= progress_text(char = '*') )
# write( names( init_rates ), file=paste( error_filehandle, "_ks.csv", sep="" ), ncolumns=length(init_rates), sep=",", append=F )
# sims_posts <- llply( .data = posts_list, .fun = function( params ){
    # init_rates[ names( params[ 1:length(init_rates) ] ) ] <- params[1:length(init_rates)]
    # init_rates <- unlist( init_rates )
    # print( init_rates )
    # # ode( y=y0, times=timesteps, func=model_eq, parms=init_rates, rtol=1e-15, method="rk4" )
    # tryCatch( ode( y=y0, times=timesteps, func=model_eq, parms=init_rates, rtol=1e-15, maxsteps=20000, method="lsoda" ),
        # warning = function(w){ print(w); write( init_rates, file=paste( error_filehandle, "_ks.csv", sep="" ), ncolumns=length(init_rates), sep=",", append=T ) } )
# }, .progress= progress_text(char = '*') )
# dead_sims <- which( sapply( sims_posts, is.null ) )
# posts_list[ dead_sims ] <- NULL

# # # Check for instances of sims where bacterial load is high for sustained period of time

# #  posts_z_list[[1]][ -match( c( "position", "combo_score" ), names( posts_z_list[[1]] ) ) ]        # remove position, combo_score
# # run_nums_2 <- lapply( posts_list, function(x){ x["position"] } )
# run_nums_2 <- lapply( posts_list, function(x){ x["position"] } )
# sims_posts <- Filter( Negate( is.null ), sims_posts )
# # run_nums_2_less <- run_nums_2[ -dead_sims ]
# # sims_posts_less <- sims_posts[ -dead_sims ]
# # sims_accept_2 <- lapply( as.list(1:length(run_nums_2_less)), function(x){ data.frame( run_num=run_nums_2_less[[x]]$position, sims_posts_less[[x]] ) } )
# # sims_accept_2 <- lapply( as.list(1:length(run_nums_2)), function(x){ data.frame( run_num=run_nums_2[[x]]$position, sims_posts[[x]] ) } )
# sims_accept_2 <- lapply( as.list(1:length(run_nums_2)), function(x){ data.frame( run_num=run_nums_2[[x]]$position, sims_posts[[x]] ) } )
# sims_accept_2 <- rbindlist( sims_accept_2 )
# sims_accept_2$log_mtb_act_gra_t <- log10( sims_accept_2$mtb_act_gra_t + 1 )
# sims_accept_2$t_cell_perc       <- 100 * sims_accept_2$imm_ant_t / ( sims_accept_2$imm_nai_t + sims_accept_2$imm_pro_t + sims_accept_2$imm_ant_t )

# threshold_gran <- 0.90
# earliest_time  <- 0        # actually, earliest time when active can occur
# weardown_dur   <- 0        # minimum duration above-threshold numbers must be maintained

# ####################

# # candidate threshold values -- look for the one that optimizes fit to the household contact data

# kamat_data <- read.csv( "kamat_data.csv" )
# kamat_data$cumul <- cumsum( kamat_data$Total )
# cumul_kamat_data <- data.frame(
    # time=c( 0, unlist( lapply( kamat_data$Time[-1], function(x){ rep(x,2) } ) )-0.5, tail( kamat_data$Time,1 )+0.5 ),
    # incid=c( unlist( lapply( rev(rev(kamat_data$cumul)), function(x){ rep(x,2) } ) ) ) )

# threshold_gran_list <- c( 0.80, 0.85, 0.90, 0.95, 0.99 )
# earliest_time_list <- c( 0, 15, 30, 60, 90, 120, 150, 180 )
# # weardown_dur_list <- c( 0, 15, 30, 60, 90, 120, 150, 180 )

# threshold_criteria <- expand.grid( threshold_gran_list, earliest_time_list )
# colnames( threshold_criteria ) <- c( "threshold_gran", "earliest_time" )
# # threshold_criteria <- expand.grid( threshold_gran_list, weardown_dur_list )
# # colnames( threshold_criteria ) <- c( "threshold_gran", "weardown_dur" )

# scores <- c()
# # k <- 1
# for ( k in 1:nrow(threshold_criteria) ) {
    
    # print( threshold_criteria[k,] )
    
    # threshold_gran_temp <- threshold_criteria$threshold_gran[k]
    # earliest_time_temp <- threshold_criteria$earliest_time[k]
    # weardown_dur_temp <- weardown_dur
    
    # sims_accept_2$above_threshold <- ( sims_accept_2$mtb_qui_gra_t + sims_accept_2$mtb_act_gra_t ) > ( threshold_gran_temp * mtb_max_gra )    # mark t when above-threshold Mtb
    # sims_accept_2$above_threshold[ which( is.na( sims_accept_2$above_threshold ) ) ] <- FALSE
    
    # active_run_pos <- c()
    # active_timepoint <- data.frame( time=c(), log_mtb_act_gra_t=c() )
    # for ( iter_2 in 1:length(unique(sims_accept_2$run_num)) ) {
        
        # if( sum( subset( sims_accept_2, run_num == (unique(sims_accept_2$run_num))[iter_2] )$above_threshold ) > 0 ) {   # if above-threshold Mtb numbers in run...
                
                # indiv_run <- subset( sims_accept_2, run_num == (unique(sims_accept_2$run_num))[iter_2] )
                
                # indiv_run_rle <- rle( indiv_run$above_threshold )
                # indiv_run_rle_df <- data.frame(
                    # start_pos = cumsum( indiv_run_rle$lengths ) - indiv_run_rle$lengths + 1,
                    # end_pos = cumsum( indiv_run_rle$lengths ),
                    # duration = indiv_run_rle$lengths,
                    # value = indiv_run_rle$values )
                
                # # check timing is after earliest time
                # meets_timing_criteria <- indiv_run$time[ indiv_run_rle_df$start_pos ] > earliest_time_temp
                
                # # check duration is longer than minimum weardown duration
                # meets_duration_criteria <- ( indiv_run$time[ indiv_run_rle_df$end_pos ] - indiv_run$time[ indiv_run_rle_df$start_pos ] ) > weardown_dur_temp
                
                # # meets all criteria (including above-threshold)?
                # indiv_run_rle_df$meets_all_criteria <- ( meets_timing_criteria * meets_duration_criteria * indiv_run_rle_df$value ) == 1
                
                # if( sum( indiv_run_rle_df$meets_all_criteria ) > 0 ) {
                    # active_run_pos <- c( active_run_pos, (unique(sims_accept_2$run_num))[iter_2] )
                    # active_timepoint <- rbind( active_timepoint, c( indiv_run[ min( indiv_run_rle_df$start_pos[ indiv_run_rle_df$meets_all_criteria ] ), c( "run_num", "time", "log_mtb_act_gra_t", "t_cell_perc" ) ] ) )
                # }
        # }
    # }
    
    # if ( length( active_run_pos ) > 0 ) {
        # incid_rle <- rle( active_timepoint[ order( active_timepoint$time ), ]$time )
        # incid_rle_df <- data.frame( time=incid_rle$values, incid=(cumsum(incid_rle$lengths))/length(posts_list) )
        
        # ann_times <- seq( from=start_time, to=end_time, by=365 )[-1]
        # ann_incid <- unlist( lapply( ann_times, function(x){ incid_rle_df$incid[ max( which( incid_rle_df$time <= x ) ) ] } ) )
        # ann_incid <- c( ann_incid[1], diff( ann_incid ) )       # http://stackoverflow.com/questions/21418287/functional-way-to-reverse-cumulative-sum
        # ann_incid_df <- data.frame( time=ann_times, incid=ann_incid )
        # ann_incid_df$time <- ann_incid_df$time/365
        # cumul_incid_rle_df <- data.frame( time=c( 0, unlist( lapply( incid_rle_df$time, function(x){ rep(x,2) } ) ), end_time ), incid=c( 0, 0, unlist( lapply( rev(rev(incid_rle_df$incid)), function(x){ rep(x,2) } ) ) ) )
        # cumul_incid_rle_df$time <- cumul_incid_rle_df$time/365
        
        # data_pts <- kamat_data[ kamat_data$Time %in% c( 1, 2, 3, 4, 5 ), "Total" ]
        # model_pts <- ann_incid_df[ ann_incid_df$time %in% c( 1, 2, 3, 4, 5 ), "incid" ]
        # scores <- c( scores, sum( unlist( lapply( ( data_pts - model_pts ), function(x){ x**2 } ) ) ) )         # sum of squares
    # } else {
        # scores <- c( scores, NA )
    # }
# }
# threshold_criteria$scores <- scores

# threshold_criteria_2dtable <- cast( threshold_criteria, earliest_time~threshold_gran, value="scores")
# # threshold_criteria_2dtable <- cast( threshold_criteria, weardown_dur~threshold_gran, value="scores")
# write.csv( threshold_criteria_2dtable, "output_fitting_hists2_threshold.csv" )

# # threshold_gran <- 0.95
# threshold_gran <- threshold_criteria$threshold_gran[ which.min( threshold_criteria$scores ) ]
# # earliest_time  <- 90        # actually, earliest time when active can occur
# earliest_time  <- threshold_criteria$earliest_time[ which.min( threshold_criteria$scores ) ]        # actually, earliest time when active can occur
# weardown_dur   <- weardown_dur        # minimum duration above-threshold numbers must be maintained
# # weardown_dur  <- threshold_criteria$weardown_dur[ which.min( threshold_criteria$scores ) ]        # actually, earliest time when active can
# print( c( threshold_gran, earliest_time ) )
# # print( c( threshold_gran, weardown_dur ) )

# #################################################################

# sims_accept_2$above_threshold <- ( sims_accept_2$mtb_qui_gra_t + sims_accept_2$mtb_act_gra_t ) > ( threshold_gran * mtb_max_gra )    # mark t when above-threshold Mtb
# sims_accept_2$above_threshold[ which( is.na( sims_accept_2$above_threshold ) ) ] <- FALSE

# active_run_pos <- c()
# active_timepoint <- data.frame( time=c(), log_mtb_act_gra_t=c() )
# for ( iter_2 in 1:length(unique(sims_accept_2$run_num)) ) {
    
    # if( sum( subset( sims_accept_2, run_num == (unique(sims_accept_2$run_num))[iter_2] )$above_threshold ) > 0 ) {   # if above-threshold Mtb numbers in run...
            
            # indiv_run <- subset( sims_accept_2, run_num == (unique(sims_accept_2$run_num))[iter_2] )
            
            # indiv_run_rle <- rle( indiv_run$above_threshold )
            # indiv_run_rle_df <- data.frame(
                # start_pos = cumsum( indiv_run_rle$lengths ) - indiv_run_rle$lengths + 1,
                # end_pos = cumsum( indiv_run_rle$lengths ),
                # duration = indiv_run_rle$lengths,
                # value = indiv_run_rle$values )
            
            # # check timing is after earliest time
            # meets_timing_criteria <- indiv_run$time[ indiv_run_rle_df$start_pos ] > earliest_time
            
            # # check duration is longer than minimum weardown duration
            # meets_duration_criteria <- ( indiv_run$time[ indiv_run_rle_df$end_pos ] - indiv_run$time[ indiv_run_rle_df$start_pos ] ) > weardown_dur
            
            # # meets all criteria (including above-threshold)?
            # indiv_run_rle_df$meets_all_criteria <- ( meets_timing_criteria * meets_duration_criteria * indiv_run_rle_df$value ) == 1
            
            # if( sum( indiv_run_rle_df$meets_all_criteria ) > 0 ) {      # if at least one instance of meeting all criteria
                # active_run_pos <- c( active_run_pos, ( unique( sims_accept_2$run_num ) )[iter_2] )
                # active_timepoint <- rbind( active_timepoint, c( indiv_run[ min( indiv_run_rle_df$start_pos[ indiv_run_rle_df$meets_all_criteria ] ), c( "run_num", "time", "log_mtb_act_gra_t", "t_cell_perc" ) ] ) )       # take the first
            # }
    # }
# }

# active_timepoint$time_rescaled <- round( active_timepoint$time * ( 1000 / max(sims_accept_2$time) ), 0 )
# active_timepoint$color <- rev( colorRampPalette( brewer.pal(9,"Reds") )( 1200 )[ 21:1020 ] )[ active_timepoint$time_rescaled ]

# ### KS #############################################

# p1_posts2 <- ggplot() +
    # theme_bw() +
    # theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
    # geom_line( data=sims_accept_2, aes( x=time, y=log_mtb_act_gra_t, group=run_num ), colour="gray20", alpha=0.5 ) +
    # geom_line( data=subset( sims_accept_2, run_num == accept_pos[ which.min( combo_score[ accept_pos ] ) ] ), aes( x=time, y=log_mtb_act_gra_t ), colour="blue", size=1.5, alpha=0.5 ) +
    # geom_line( data=sims_accept_2[ sims_accept_2$run_num %in% active_run_pos, ], aes( x=time, y=log_mtb_act_gra_t, colour=factor(run_num) ), size=1.5, alpha=0.5 ) +
    # geom_point( data=active_timepoint, aes( x=time, y=log_mtb_act_gra_t ), pch=13, size=3 ) +
    # # geom_point( data=active_timepoint, aes( x=time, y=log_mtb_act_gra_t, colour=factor(run_num) ), pch=13, size=3 ) +
    # xlab( "Time (days)" ) +
    # ylab( "log10 CFU" ) +
    # scale_color_manual( values=active_timepoint$color ) +
    # theme( legend.position="none" ) +
    # ylim( c( 0, log10( mtb_max_gra + 1 ) ) )

# p1_posts3 <- ggplot() +
    # theme_bw() +
    # theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
    # geom_line( data=sims_accept_2, aes( x=time, y=t_cell_perc, group=run_num ), colour="gray20", alpha=0.5 ) +
    # geom_line( data=subset( sims_accept_2, run_num == accept_pos[ which.min( combo_score[ accept_pos ] ) ] ), aes( x=time, y=t_cell_perc ), colour="blue", size=1.5, alpha=0.5 ) +
    # geom_line( data=sims_accept_2[ sims_accept_2$run_num %in% active_run_pos, ], aes( x=time, y=t_cell_perc, colour=factor(run_num) ), size=1.5, alpha=0.5 ) +
    # geom_point( data=active_timepoint, aes( x=time, y=t_cell_perc ), pch=13, size=3 ) +
    # # geom_point( data=active_timepoint, aes( x=time, y=log_mtb_act_gra_t, colour=factor(run_num) ), pch=13, size=3 ) +
    # xlab( "Time (days)" ) +
    # ylab( "Treg %age of CD4" ) +
    # ylim( c( 0, 100 ) ) +
    # scale_color_manual( values=active_timepoint$color ) +
    # theme( legend.position="none" )

# incid_rle <- rle( active_timepoint[ order( active_timepoint$time ), ]$time )
# incid_rle_df <- data.frame( time=incid_rle$values, incid=(cumsum(incid_rle$lengths))/length(posts_list) )

# ann_times <- seq( from=start_time, to=end_time, by=365 )[-1]
# ann_incid <- unlist( lapply( ann_times, function(x){ incid_rle_df$incid[ max( which( incid_rle_df$time <= x ) ) ] } ) )
# ann_incid <- c( ann_incid[1], diff( ann_incid ) )       # http://stackoverflow.com/questions/21418287/functional-way-to-reverse-cumulative-sum
# ann_incid_df <- data.frame( time=ann_times, incid=ann_incid )
# cumul_incid_rle_df <- data.frame( time=c( 0, unlist( lapply( incid_rle_df$time, function(x){ rep(x,2) } ) ), end_time ), incid=c( 0, 0, unlist( lapply( rev(rev(incid_rle_df$incid)), function(x){ rep(x,2) } ) ) ) )

# # incid_plot <- ggplot() +
    # # theme_bw() +
    # # theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
    # # # geom_point( data=incid_alt_rle_df, aes( x=time, y=incid ), pch=13, size=3 ) +
    # # geom_line( data=cumul_incid_rle_df, aes( x=time, y=incid ), color="blue" ) +
    # # geom_point( data=ann_incid_df, aes( x=time, y=incid ), pch=16, size=3, color="blue" ) +
    # # xlab( "Time (days)" ) +
    # # ylab( "Incidence" ) +
    # # xlim( c( 0, end_time ) )

# kamat_data <- read.csv( "kamat_data.csv" )
# kamat_data$cumul <- cumsum( kamat_data$Total )
# cumul_kamat_data <- data.frame(
    # time=c( 0, unlist( lapply( kamat_data$Time[-1], function(x){ rep(x,2) } ) )-0.5, tail( kamat_data$Time,1 )+0.5 ),
    # incid=c( unlist( lapply( rev(rev(kamat_data$cumul)), function(x){ rep(x,2) } ) ) ) )

# incid_plot <- ggplot( ) +
    # theme_bw() +
    # theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
    # # geom_point( data=incid_alt_rle_df, aes( x=time, y=incid ), pch=13, size=4 ) +
    # geom_line( data=cumul_incid_rle_df, aes( x=time/365, y=incid*100, color="Model" ), size=1.5, alpha=0.8 ) +
    # geom_point( data=ann_incid_df, aes( x=time/365, y=incid*100, color="Model"), pch=15, size=3.5, alpha=0.6 ) +
    # geom_point( data=ann_incid_df, aes( x=time/365, y=incid*100 ), pch=0, size=3.5, color="black", alpha=0.6 ) +
    # xlab( "Time (years)" ) +
    # ylab( "Incidence (% of population)" ) +
    # # xlim( c( 0, end_time/365 ) ) +
    # # geom_line( data=kamat_data, aes( x=Time, y=cumul ), color="red" ) +
    # geom_line( data=cumul_kamat_data, aes( x=time, y=incid*100, color="Kamat" ), size=1.5, alpha=0.8 ) +
    # geom_point( data=kamat_data, aes( x=Time, y=Total*100, color="Kamat" ), pch=19, size=4, alpha=0.6 ) +
    # geom_point( data=kamat_data, aes( x=Time, y=Total*100 ), pch=1, size=4, color="black", alpha=0.6 ) +
    # scale_x_continuous( breaks = c( 0, 2, 4, 6, 8, 10 ) ) +
    # # scale_y_continuous( breaks = c( 0, 2, 4, 6, 8, 10, 12 ) ) +
    # scale_y_continuous( breaks = seq( 0, round( 100*max( cumul_incid_rle_df$incid, cumul_kamat_data$y ) ), 2 ) ) +
    # scale_color_manual( name="Source", values=c( Model="blue", Kamat="red" ) ) +
    # labs( caption = paste( "threshold=", threshold_gran, ", earliest time=", earliest_time, sep="" ) )

# pdf( height=5, width=8, "output_fitting_hists2.pdf" )
# print( p1_posts2 )
# print( p1_posts3 )
# print( incid_plot )
# dev.off()

# write.csv( sims_accept_2, "sims_accept.csv", row.names=F )
# write.csv( active_timepoint, "active_timepoint.csv", row.names=F )


###################################################
### Z #############################################
###################################################

posts_z_list <- vector( mode="list", length=nrow( posts_z ) )
for ( i in 1:length( posts_z_list ) ) {
    # print( i )
    posts_z_list[[ i ]] <- posts_z[ i, ]
}

# sims_posts_z <- llply( .data = posts_z_list, .fun = function( lhs_params ){
    # init_rates[ ( names( lhs_params )[1:length(init_rates)] ) ] <- lhs_params[1:length(init_rates)]
    # init_rates <- unlist( init_rates )
    # print( init_rates )
    # ode( y=y0, times=timesteps, func=model_eq, parms=init_rates, rtol=1e-15, maxsteps=20000, method="lsoda" )
# }, .progress= progress_text(char = '*') )
write( names( init_rates ), file=paste( error_filehandle, "_z.csv", sep="" ), ncolumns=length(init_rates), sep=",", append=F )
sims_posts_z <- llply( .data = posts_z_list, .fun = function( params ){
    init_rates[ names( params[ 1:length(init_rates) ] ) ] <- params[1:length(init_rates)]
    init_rates <- unlist( init_rates )
    print( init_rates )
    # ode( y=y0, times=timesteps, func=model_eq, parms=init_rates, rtol=1e-15, method="rk4" )
    # tryCatch( ode( y=y0, times=timesteps, func=model_eq, parms=init_rates, rtol=1e-15, maxsteps=20000, method="lsoda" ),
    tryCatch( ode( y=y0, times=timesteps, func=model_eq, parms=init_rates, method="euler", hini=1 ),
        warning = function(w){ print(w); write( init_rates, file=paste( error_filehandle, "_z.csv", sep="" ), ncolumns=length(init_rates), sep=",", append=T ) } )
}, .progress= progress_text(char = '*') )
dead_sims_z <- which( sapply( sims_posts_z, is.null ) )
posts_z_list[ dead_sims_z ] <- NULL

# # Check for instances of sims where bacterial load is high for sustained period of time

#  posts_z_list[[1]][ -match( c( "position", "combo_score" ), names( posts_z_list[[1]] ) ) ]        # remove position, combo_score
# run_nums_z2 <- lapply( posts_z_list, function(x){ x["position"] } )
run_nums_z2 <- lapply( posts_z_list, function(x){ x["position"] } )
sims_posts_z <- Filter( Negate( is.null ), sims_posts_z )
# sims_accept_z2 <- lapply( as.list(1:length(run_nums_z2_less)), function(x){ data.frame( run_num=run_nums_z2_less[[x]]$position, sims_posts_z_less[[x]] ) } )
# sims_accept_z2 <- lapply( as.list(1:length(run_nums_z2)), function(x){ if( !is.null( sims_posts_z[[x]] ) ) { data.frame( run_num=run_nums_z2[[x]]$position, sims_posts_z[[x]] ) } } )
sims_accept_z2 <- lapply( as.list(1:length(run_nums_z2)), function(x){ data.frame( run_num=run_nums_z2[[x]]$position, sims_posts_z[[x]] ) } )
sims_accept_z2 <- rbindlist( sims_accept_z2 )
sims_accept_z2$log_mtb_act_gra_t <- log10( sims_accept_z2$mtb_act_gra_t + 1 )
sims_accept_z2$t_cell_perc       <- 100 * sims_accept_z2$imm_ant_t / ( sims_accept_z2$imm_nai_t + sims_accept_z2$imm_pro_t + sims_accept_z2$imm_ant_t )

threshold_gran_z <- 0.90
earliest_time_z  <- 0        # actually, earliest time when active can occur
weardown_dur_z   <- 0        # minimum duration above-threshold numbers must be maintained

####################

# candidate threshold values -- look for the one that optimizes fit to the household contact data

kamat_data <- read.csv( "kamat_data.csv" )
kamat_data$cumul <- cumsum( kamat_data$Total )
cumul_kamat_data <- data.frame(
    time=c( 0, unlist( lapply( kamat_data$Time[-1], function(x){ rep(x,2) } ) )-0.5, tail( kamat_data$Time,1 )+0.5 ),
    incid=c( unlist( lapply( rev(rev(kamat_data$cumul)), function(x){ rep(x,2) } ) ) ) )

threshold_gran_list <- c( 0.80, 0.85, 0.90, 0.95, 0.99 )
earliest_time_list <- c( 0, 15, 30, 60, 90, 120, 150, 180 )
# expand.grid( c(1,2,3), c("A","B") )
threshold_criteria <- expand.grid( threshold_gran_list, earliest_time_list )
colnames( threshold_criteria ) <- c( "threshold_gran", "earliest_time" )

scores <- c()
# k <- 1
for ( k in 1:nrow(threshold_criteria) ) {
    
    print( threshold_criteria[k,] )
    
    threshold_gran_temp <- threshold_criteria$threshold_gran[k]
    earliest_time_temp <- threshold_criteria$earliest_time[k]
    weardown_dur_temp <- weardown_dur_z
    
    sims_accept_z2$above_threshold <- ( sims_accept_z2$mtb_qui_gra_t + sims_accept_z2$mtb_act_gra_t ) > ( threshold_gran_temp * mtb_max_gra )    # mark t when above-threshold Mtb
    sims_accept_z2$above_threshold[ which( is.na( sims_accept_z2$above_threshold ) ) ] <- FALSE
    
    active_run_pos_z <- c()
    active_timepoint_z <- data.frame( time=c(), log_mtb_act_gra_t=c() )
    for ( iter_z2 in 1:length(unique(sims_accept_z2$run_num)) ) {
        
        if( sum( subset( sims_accept_z2, run_num == (unique(sims_accept_z2$run_num))[iter_z2] )$above_threshold ) > 0 ) {   # if above-threshold Mtb numbers in run...
                
                indiv_run <- subset( sims_accept_z2, run_num == (unique(sims_accept_z2$run_num))[iter_z2] )
                
                indiv_run_rle <- rle( indiv_run$above_threshold )
                indiv_run_rle_df <- data.frame(
                    start_pos = cumsum( indiv_run_rle$lengths ) - indiv_run_rle$lengths + 1,
                    end_pos = cumsum( indiv_run_rle$lengths ),
                    duration = indiv_run_rle$lengths,
                    value = indiv_run_rle$values )
                
                # check timing is after earliest time
                meets_timing_criteria <- indiv_run$time[ indiv_run_rle_df$start_pos ] > earliest_time_temp
                
                # check duration is longer than minimum weardown duration
                meets_duration_criteria <- ( indiv_run$time[ indiv_run_rle_df$end_pos ] - indiv_run$time[ indiv_run_rle_df$start_pos ] ) > weardown_dur_z
                
                # meets all criteria (including above-threshold)?
                indiv_run_rle_df$meets_all_criteria <- ( meets_timing_criteria * meets_duration_criteria * indiv_run_rle_df$value ) == 1
                
                if( sum( indiv_run_rle_df$meets_all_criteria ) > 0 ) {
                    active_run_pos_z <- c( active_run_pos_z, (unique(sims_accept_z2$run_num))[iter_z2] )
                    active_timepoint_z <- rbind( active_timepoint_z, c( indiv_run[ min( indiv_run_rle_df$start_pos[ indiv_run_rle_df$meets_all_criteria ] ), c( "run_num", "time", "log_mtb_act_gra_t", "t_cell_perc" ) ] ) )
                }
        }
    }
    
    if ( length( active_run_pos_z ) > 0 ) {
        incid_z_rle <- rle( active_timepoint_z[ order( active_timepoint_z$time ), ]$time )
        incid_z_rle_df <- data.frame( time=incid_z_rle$values, incid=(cumsum(incid_z_rle$lengths))/length(posts_z_list) )
        
        ann_times <- seq( from=start_time, to=end_time, by=365 )[-1]
        ann_incid <- unlist( lapply( ann_times, function(x){ incid_z_rle_df$incid[ max( which( incid_z_rle_df$time <= x ) ) ] } ) )
        ann_incid <- c( ann_incid[1], diff( ann_incid ) )       # http://stackoverflow.com/questions/21418287/functional-way-to-reverse-cumulative-sum
        ann_incid_df <- data.frame( time=ann_times, incid=ann_incid )
        ann_incid_df$time <- ann_incid_df$time/365
        # cumul_incid_z_rle_df <- data.frame( time=c( 0, unlist( lapply( incid_z_rle_df$time, function(x){ rep(x,2) } ) ), end_time ), incid=c( 0, 0, unlist( lapply( rev(rev(incid_z_rle_df$incid)), function(x){ rep(x,2) } ) ) ) )
        # cumul_incid_z_rle_df$time <- cumul_incid_z_rle_df$time/365
        
        data_pts <- kamat_data[ kamat_data$Time %in% c( 1, 2, 3, 4, 5 ), "Total" ]
        model_pts <- ann_incid_df[ ann_incid_df$time %in% c( 1, 2, 3, 4, 5 ), "incid" ]
        scores <- c( scores, sum( unlist( lapply( ( data_pts - model_pts ), function(x){ x**2 } ) ) ) )         # sum of squares
    } else {
        scores <- c( scores, NA )
    }
}
threshold_criteria$scores <- scores
threshold_criteria_2dtable <- cast( threshold_criteria, earliest_time~threshold_gran, value="scores")
write.csv( threshold_criteria_2dtable, "output_fitting_z_hists2_threshold.csv" )

threshold_gran_z <- threshold_criteria$threshold_gran[ which.min( threshold_criteria$scores ) ]
earliest_time_z  <- threshold_criteria$earliest_time[ which.min( threshold_criteria$scores ) ]        # actually, earliest time when active can occur
weardown_dur_z   <- weardown_dur_z        # minimum duration above-threshold numbers must be maintained
print( c( threshold_gran_z, earliest_time_z ) )

######################################

sims_accept_z2$above_threshold <- ( sims_accept_z2$mtb_qui_gra_t + sims_accept_z2$mtb_act_gra_t ) > ( threshold_gran_z * mtb_max_gra )    # mark t when above-threshold Mtb
sims_accept_z2$above_threshold[ which( is.na( sims_accept_z2$above_threshold ) ) ] <- FALSE

active_run_pos_z <- c()
active_timepoint_z <- data.frame( time=c(), log_mtb_act_gra_t=c() )
for ( iter_z2 in 1:length(unique(sims_accept_z2$run_num)) ) {
    
    if( sum( subset( sims_accept_z2, run_num == (unique(sims_accept_z2$run_num))[iter_z2] )$above_threshold ) > 0 ) {   # if above-threshold Mtb numbers in run...
            
            indiv_run <- subset( sims_accept_z2, run_num == (unique(sims_accept_z2$run_num))[iter_z2] )
            
            indiv_run_rle <- rle( indiv_run$above_threshold )
            indiv_run_rle_df <- data.frame(
                start_pos = cumsum( indiv_run_rle$lengths ) - indiv_run_rle$lengths + 1,
                end_pos = cumsum( indiv_run_rle$lengths ),
                duration = indiv_run_rle$lengths,
                value = indiv_run_rle$values )
            
            # check timing is after earliest time
            meets_timing_criteria <- indiv_run$time[ indiv_run_rle_df$start_pos ] > earliest_time_z
            
            # check duration is longer than minimum weardown duration
            meets_duration_criteria <- ( indiv_run$time[ indiv_run_rle_df$end_pos ] - indiv_run$time[ indiv_run_rle_df$start_pos ] ) >= weardown_dur_z
            
            # meets all criteria (including above-threshold)?
            indiv_run_rle_df$meets_all_criteria <- ( meets_timing_criteria * meets_duration_criteria * indiv_run_rle_df$value ) == 1
            
            if( sum( indiv_run_rle_df$meets_all_criteria ) > 0 ) {
                active_run_pos_z <- c( active_run_pos_z, ( unique( sims_accept_z2$run_num ) )[iter_z2] )
                active_timepoint_z <- rbind( active_timepoint_z, c( indiv_run[ min( indiv_run_rle_df$start_pos[ indiv_run_rle_df$meets_all_criteria ] ), c( "run_num", "time", "log_mtb_act_gra_t", "t_cell_perc" ) ] ) )
            }
    }
}

active_timepoint_z$time_rescaled <- round( active_timepoint_z$time * ( 1000 / max(sims_accept_z2$time) ), 0 )
active_timepoint_z$color <- rev( colorRampPalette( brewer.pal(9,"Reds") )( 1200 )[ 21:1020 ] )[ active_timepoint_z$time_rescaled ]

### Z #############################################

p1_z_posts2 <- ggplot() +
    theme_bw() +
    theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
    geom_line( data=sims_accept_z2, aes( x=time, y=log_mtb_act_gra_t, group=run_num ), colour="gray20", alpha=0.5 ) +
    geom_line( data=subset( sims_accept_z2, run_num == accept_pos_z[ which.min( combo_score_z[ accept_pos_z ] ) ] ), aes( x=time, y=log_mtb_act_gra_t ), colour="blue", size=1.5, alpha=0.5 ) +
    geom_line( data=sims_accept_z2[ sims_accept_z2$run_num %in% active_run_pos_z, ], aes( x=time, y=log_mtb_act_gra_t, colour=factor(run_num) ), size=1.5, alpha=0.5 ) +
    geom_point( data=active_timepoint_z, aes( x=time, y=log_mtb_act_gra_t ), pch=13, size=3 ) +
    # geom_point( data=active_timepoint, aes( x=time, y=log_mtb_act_gra_t, colour=factor(run_num) ), pch=13, size=3 ) +
    xlab( "Time (days)" ) +
    ylab( "log10 CFU" ) +
    scale_color_manual( values=active_timepoint_z$color ) +
    theme( legend.position="none" ) +
    ylim( c( 0, log10( mtb_max_gra + 1 ) ) )

p1_z_posts3 <- ggplot() +
    theme_bw() +
    theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
    geom_line( data=sims_accept_z2, aes( x=time, y=t_cell_perc, group=run_num ), colour="gray20", alpha=0.5 ) +
    geom_line( data=subset( sims_accept_z2, run_num == accept_pos_z[ which.min( combo_score_z[ accept_pos_z ] ) ] ), aes( x=time, y=t_cell_perc ), colour="blue", size=1.5, alpha=0.5 ) +
    geom_line( data=sims_accept_z2[ sims_accept_z2$run_num %in% active_run_pos_z, ], aes( x=time, y=t_cell_perc, colour=factor(run_num) ), size=1.5, alpha=0.5 ) +
    geom_point( data=active_timepoint_z, aes( x=time, y=t_cell_perc ), pch=13, size=3 ) +
    # geom_point( data=active_timepoint, aes( x=time, y=log_mtb_act_gra_t, colour=factor(run_num) ), pch=13, size=3 ) +
    xlab( "Time (days)" ) +
    ylab( "Treg %age of CD4" ) +
    ylim( c( 0, 100 ) ) +
    scale_color_manual( values=active_timepoint_z$color ) +
    theme( legend.position="none" )

incid_z_rle <- rle( active_timepoint_z[ order( active_timepoint_z$time ), ]$time )
incid_z_rle_df <- data.frame( time=incid_z_rle$values, incid=(cumsum(incid_z_rle$lengths))/length(posts_z_list) )

# incid_z_plot <- ggplot() +
    # theme_bw() +
    # theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
    # geom_point( data=incid_z_rle_df, aes( x=time, y=incid ), pch=13, size=3 ) +
    # xlab( "Time (days)" ) +
    # ylab( "Cumulative incidence" ) +
    # xlim( c( 0, max( sims_accept_z2$time ) ) )

ann_times <- seq( from=start_time, to=end_time, by=365 )[-1]
ann_incid <- unlist( lapply( ann_times, function(x){ incid_z_rle_df$incid[ max( which( incid_z_rle_df$time <= x ) ) ] } ) )
ann_incid <- c( ann_incid[1], diff( ann_incid ) )       # http://stackoverflow.com/questions/21418287/functional-way-to-reverse-cumulative-sum
ann_incid_df <- data.frame( time=ann_times, incid=ann_incid )
# cum_incid_z_rle_df <- data.frame( time=c( 0, unlist( lapply( incid_z_rle_df$time, function(x){ rep(x,2) } ) ) ), incid=c( 0, 0, unlist(lapply( rev(rev(incid_z_rle_df$incid)[-1]), function(x){ rep(x,2) } )), rev(incid_z_rle_df$incid)[1] ) )
cumul_incid_z_rle_df <- data.frame( time=c( 0, unlist( lapply( incid_z_rle_df$time, function(x){ rep(x,2) } ) ), end_time ), incid=c( 0, 0, unlist( lapply( rev(rev(incid_z_rle_df$incid)), function(x){ rep(x,2) } ) ) ) )

# incid_z_plot <- ggplot() +
    # theme_bw() +
    # theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
    # # geom_point( data=incid_alt_rle_df, aes( x=time, y=incid ), pch=13, size=3 ) +
    # geom_line( data=cumul_incid_z_rle_df, aes( x=time, y=incid ), color="blue" ) +
    # geom_point( data=ann_incid_df, aes( x=time, y=incid ), pch=16, size=3, color="blue" ) +
    # xlab( "Time (days)" ) +
    # ylab( "Incidence" ) +
    # xlim( c( 0, end_time ) )

kamat_data <- read.csv( "kamat_data.csv" )
kamat_data$cumul <- cumsum( kamat_data$Total )
cumul_kamat_data <- data.frame(
    time=c( 0, unlist( lapply( kamat_data$Time[-1], function(x){ rep(x,2) } ) )-0.5, tail( kamat_data$Time,1 )+0.5 ),
    incid=c( unlist( lapply( rev(rev(kamat_data$cumul)), function(x){ rep(x,2) } ) ) ) )

incid_z_plot <- ggplot() +
    theme_bw() +
    theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
    # geom_line( data=cum_incid_z_rle_df, aes( x=time, y=incid ), color="blue" ) +
    geom_line( data=cumul_incid_z_rle_df, aes( x=time/365, y=incid*100, color="Model" ), size=1.5, alpha=0.8 ) +
    geom_point( data=ann_incid_df, aes( x=time/365, y=incid*100, color="Model"), pch=15, size=3.5, alpha=0.6 ) +
    geom_point( data=ann_incid_df, aes( x=time/365, y=incid*100 ), pch=0, size=3.5, color="black", alpha=0.6 ) +
    # geom_point( data=ann_incid_df, aes( x=time, y=incid ), pch=16, size=3, color="blue" ) +
    xlab( "Time (years)" ) +
    ylab( "Incidence (% of population)" ) +
    # xlim( c( 0, end_time ) )
    geom_line( data=cumul_kamat_data, aes( x=time, y=incid*100, color="Kamat" ), size=1.5, alpha=0.8 ) +
    geom_point( data=kamat_data, aes( x=Time, y=Total*100, color="Kamat" ), pch=19, size=4, alpha=0.6 ) +
    geom_point( data=kamat_data, aes( x=Time, y=Total*100 ), pch=1, size=4, color="black", alpha=0.6 ) +
    scale_x_continuous( breaks = c( 0, 2, 4, 6, 8, 10 ) ) +
    # scale_y_continuous( breaks = c( 0, 2, 4, 6, 8, 10, 12 ) ) +
    scale_y_continuous( breaks = seq( 0, round( 100*max( cumul_incid_z_rle_df$incid, cumul_kamat_data$y ) ), 2 ) ) +
    scale_color_manual( name="Source", values=c( Model="blue", Kamat="red" ) ) +
    labs( caption = paste( "threshold=", threshold_gran_z, ", earliest time=", earliest_time_z, sep="" ) )

pdf( height=5, width=8, "output_fitting_z_hists2.pdf" )
print( p1_z_posts2 )
print( p1_z_posts3 )
print( incid_z_plot )
dev.off()

write.csv( sims_accept_z2, "sims_accept_z.csv", row.names=F )
write.csv( active_timepoint_z, "active_timepoint_z.csv", row.names=F )


#####################################################
### alt #############################################
#####################################################

# y0         <- c( mtb_qui_gra_t=mtb_qui_gra_init, mtb_act_gra_t=mtb_act_gra_init, imm_nai_t=imm_nai_init, imm_pro_t=imm_pro_init, imm_ant_t=imm_ant_init )

# init_rates <- c( growth_rate_log=growth_rate_log, cmi_kill_rate_log=cmi_kill_rate_log, imm_source_rate_log=imm_source_rate_log, nai2pro_diff_rate_log=nai2pro_diff_rate_log, nai2ant_diff_rate_log=nai2ant_diff_rate_log, nai_decay_rate_log=nai_decay_rate_log, pro_decay_rate_log=pro_decay_rate_log, ant_decay_rate_log=ant_decay_rate_log )

posts_alt_list <- vector( mode="list", length=nrow( posts_alt ) )
for ( i in 1:length( posts_alt_list ) ) {
    # print( i )
    posts_alt_list[[ i ]] <- posts_alt[ i, ]
}

# sims_posts_alt <- llply( .data = posts_alt_list, .fun = function( lhs_params ){
    # init_rates[ ( names( lhs_params )[1:length(init_rates)] ) ] <- lhs_params[1:length(init_rates)]
    # # init_rates[ names( lhs_params[ 1:length(init_rates) ] ) ] <- lhs_params
    # init_rates <- unlist( init_rates )
    # print( init_rates )
    # ode( y=y0, times=timesteps, func=model_eq, parms=init_rates, rtol=1e-15, maxsteps=20000, method="lsoda" )
# }, .progress= progress_text(char = '*') )
write( names( init_rates ), file=paste( error_filehandle, "_alt.csv", sep="" ), ncolumns=length(init_rates), sep=",", append=F )
sims_posts_alt <- llply( .data = posts_alt_list, .fun = function( params ){
    init_rates[ names( params[ 1:length(init_rates) ] ) ] <- params[1:length(init_rates)]
    init_rates <- unlist( init_rates )
    print( init_rates )
    # ode( y=y0, times=timesteps, func=model_eq, parms=init_rates, rtol=1e-15, method="rk4" )
    # tryCatch( ode( y=y0, times=timesteps, func=model_eq, parms=init_rates, rtol=1e-15, maxsteps=20000, method="lsoda" ),
    tryCatch( ode( y=y0, times=timesteps, func=model_eq, parms=init_rates, method="euler", hini=1 ),
        warning = function(w){ print(w); write( init_rates, file=paste( error_filehandle, "_alt.csv", sep="" ), ncolumns=length(init_rates), sep=",", append=T ) } )
}, .progress= progress_text(char = '*') )
dead_sims_alt <- which( sapply( sims_posts_alt, is.null ) )
posts_alt_list[ dead_sims_alt ] <- NULL

# # Check for instances of sims where bacterial load is high for sustained period of time

#  posts_z_list[[1]][ -match( c( "position", "combo_score" ), names( posts_z_list[[1]] ) ) ]        # remove position, combo_score
# run_nums_alt2 <- lapply( posts_alt_list, function(x){ x["position"] } )
run_nums_alt2 <- lapply( posts_alt_list, function(x){ x["position"] } )
sims_posts_alt <- Filter( Negate( is.null ), sims_posts_alt )
# sims_accept_alt2 <- lapply( as.list(1:length(run_nums_alt2_less)), function(x){ data.frame( run_num=run_nums_alt2_less[[x]]$position, sims_posts_alt_less[[x]] ) } )
# sims_accept_alt2 <- lapply( as.list(1:length(run_nums_alt2)), function(x){ data.frame( run_num=run_nums_alt2[[x]]$position, sims_posts_alt[[x]] ) } )
# sims_accept_alt2 <- lapply( as.list(1:length(run_nums_alt2)), function(x){ if( !is.null( sims_posts_alt2[[x]] ) ) { data.frame( run_num=run_nums_alt2[[x]]$position, sims_posts_alt[[x]] ) } } )
sims_accept_alt2 <- lapply( as.list(1:length(run_nums_alt2)), function(x){ data.frame( run_num=run_nums_alt2[[x]]$position, sims_posts_alt[[x]] ) } )
sims_accept_alt2 <- rbindlist( sims_accept_alt2 )
sims_accept_alt2$log_mtb_act_gra_t <- log10( sims_accept_alt2$mtb_act_gra_t + 1 )
sims_accept_alt2$t_cell_perc       <- 100 * sims_accept_alt2$imm_ant_t / ( sims_accept_alt2$imm_nai_t + sims_accept_alt2$imm_pro_t + sims_accept_alt2$imm_ant_t )

threshold_gran_alt <- 0.90
earliest_time_alt  <- 0        # actually, earliest time when active can occur
weardown_dur_alt   <- 0        # minimum duration above-threshold numbers must be maintained

####################

# candidate threshold values -- look for the one that optimizes fit to the household contact data

kamat_data <- read.csv( "kamat_data.csv" )
kamat_data$cumul <- cumsum( kamat_data$Total )
cumul_kamat_data <- data.frame(
    time=c( 0, unlist( lapply( kamat_data$Time[-1], function(x){ rep(x,2) } ) )-0.5, tail( kamat_data$Time,1 )+0.5 ),
    incid=c( unlist( lapply( rev(rev(kamat_data$cumul)), function(x){ rep(x,2) } ) ) ) )

threshold_gran_list <- c( 0.80, 0.85, 0.90, 0.95, 0.99 )
earliest_time_list <- c( 0, 15, 30, 60, 90, 120, 150, 180 )
# expand.grid( c(1,2,3), c("A","B") )
threshold_criteria <- expand.grid( threshold_gran_list, earliest_time_list )
colnames( threshold_criteria ) <- c( "threshold_gran", "earliest_time" )

scores <- c()
# k <- 1
for ( k in 1:nrow(threshold_criteria) ) {
    
    print( threshold_criteria[k,] )
    
    threshold_gran_temp <- threshold_criteria$threshold_gran[k]
    earliest_time_temp <- threshold_criteria$earliest_time[k]
    weardown_dur_temp <- weardown_dur_alt
    
    sims_accept_alt2$above_threshold <- ( sims_accept_alt2$mtb_qui_gra_t + sims_accept_alt2$mtb_act_gra_t ) > ( threshold_gran_temp * mtb_max_gra )    # mark t when above-threshold Mtb
    sims_accept_alt2$above_threshold[ which( is.na( sims_accept_alt2$above_threshold ) ) ] <- FALSE
    
    active_run_pos_alt <- c()
    active_timepoint_alt <- data.frame( time=c(), log_mtb_act_gra_t=c() )
    for ( iter_alt2 in 1:length(unique(sims_accept_alt2$run_num)) ) {
        
        if( sum( subset( sims_accept_alt2, run_num == (unique(sims_accept_alt2$run_num))[iter_alt2] )$above_threshold ) > 0 ) {   # if above-threshold Mtb numbers in run...
                
                indiv_run <- subset( sims_accept_alt2, run_num == (unique(sims_accept_alt2$run_num))[iter_alt2] )
                
                indiv_run_rle <- rle( indiv_run$above_threshold )
                indiv_run_rle_df <- data.frame(
                    start_pos = cumsum( indiv_run_rle$lengths ) - indiv_run_rle$lengths + 1,
                    end_pos = cumsum( indiv_run_rle$lengths ),
                    duration = indiv_run_rle$lengths,
                    value = indiv_run_rle$values )
                
                # check timing is after earliest time
                meets_timing_criteria <- indiv_run$time[ indiv_run_rle_df$start_pos ] > earliest_time_temp
                
                # check duration is longer than minimum weardown duration
                meets_duration_criteria <- ( indiv_run$time[ indiv_run_rle_df$end_pos ] - indiv_run$time[ indiv_run_rle_df$start_pos ] ) > weardown_dur_temp
                
                # meets all criteria (including above-threshold)?
                indiv_run_rle_df$meets_all_criteria <- ( meets_timing_criteria * meets_duration_criteria * indiv_run_rle_df$value ) == 1
                
                if( sum( indiv_run_rle_df$meets_all_criteria ) > 0 ) {
                    active_run_pos_alt <- c( active_run_pos_alt, (unique(sims_accept_alt2$run_num))[iter_alt2] )
                    active_timepoint_alt <- rbind( active_timepoint_alt, c( indiv_run[ min( indiv_run_rle_df$start_pos[ indiv_run_rle_df$meets_all_criteria ] ), c( "run_num", "time", "log_mtb_act_gra_t", "t_cell_perc" ) ] ) )
                }
        }
    }
    
    if ( length( active_run_pos_alt ) > 0 ) {
        incid_alt_rle <- rle( active_timepoint_alt[ order( active_timepoint_alt$time ), ]$time )
        incid_alt_rle_df <- data.frame( time=incid_alt_rle$values, incid=(cumsum(incid_alt_rle$lengths))/length(posts_alt_list) )
        
        ann_times <- seq( from=start_time, to=end_time, by=365 )[-1]
        ann_incid <- unlist( lapply( ann_times, function(x){ incid_alt_rle_df$incid[ max( which( incid_alt_rle_df$time <= x ) ) ] } ) )
        ann_incid <- c( ann_incid[1], diff( ann_incid ) )       # http://stackoverflow.com/questions/21418287/functional-way-to-reverse-cumulative-sum
        ann_incid_df <- data.frame( time=ann_times, incid=ann_incid )
        ann_incid_df$time <- ann_incid_df$time/365
        cumul_incid_alt_rle_df <- data.frame( time=c( 0, unlist( lapply( incid_alt_rle_df$time, function(x){ rep(x,2) } ) ), end_time ), incid=c( 0, 0, unlist( lapply( rev(rev(incid_alt_rle_df$incid)), function(x){ rep(x,2) } ) ) ) )
        cumul_incid_alt_rle_df$time <- cumul_incid_alt_rle_df$time/365
        
        data_pts <- kamat_data[ kamat_data$Time %in% c( 1, 2, 3, 4, 5 ), "Total" ]
        model_pts <- ann_incid_df[ ann_incid_df$time %in% c( 1, 2, 3, 4, 5 ), "incid" ]
        scores <- c( scores, sum( unlist( lapply( data_pts-model_pts, function(x){ x**2 } ) ) ) )         # sum of squares
    } else {
        scores <- c( scores, NA )
    }
}
threshold_criteria$scores <- scores
threshold_criteria_2dtable <- cast( threshold_criteria, earliest_time~threshold_gran, value="scores")
write.csv( threshold_criteria_2dtable, "output_fitting_alt_hists2_threshold.csv" )

# chosen threshold values below

threshold_gran_alt <- threshold_criteria$threshold_gran[ which.min( threshold_criteria$scores ) ]
earliest_time_alt  <- threshold_criteria$earliest_time[ which.min( threshold_criteria$scores ) ]        # actually, earliest time when active can occur
weardown_dur_alt   <- weardown_dur_alt        # minimum duration above-threshold numbers must be maintained
print( c( threshold_gran_alt, earliest_time_alt ) )

####################

sims_accept_alt2$above_threshold <- ( sims_accept_alt2$mtb_qui_gra_t + sims_accept_alt2$mtb_act_gra_t ) > ( threshold_gran_alt * mtb_max_gra )    # mark t when above-threshold Mtb
sims_accept_alt2$above_threshold[ which( is.na( sims_accept_alt2$above_threshold ) ) ] <- FALSE

active_run_pos_alt <- c()
active_timepoint_alt <- data.frame( time=c(), log_mtb_act_gra_t=c() )
for ( iter_alt2 in 1:length(unique(sims_accept_alt2$run_num)) ) {
    
    if( sum( subset( sims_accept_alt2, run_num == (unique(sims_accept_alt2$run_num))[iter_alt2] )$above_threshold ) > 0 ) {   # if above-threshold Mtb numbers in run...
            
            indiv_run <- subset( sims_accept_alt2, run_num == (unique(sims_accept_alt2$run_num))[iter_alt2] )
            
            indiv_run_rle <- rle( indiv_run$above_threshold )
            indiv_run_rle_df <- data.frame(
                start_pos = cumsum( indiv_run_rle$lengths ) - indiv_run_rle$lengths + 1,
                end_pos = cumsum( indiv_run_rle$lengths ),
                duration = indiv_run_rle$lengths,
                value = indiv_run_rle$values )
            
            # check timing is after earliest time
            meets_timing_criteria <- indiv_run$time[ indiv_run_rle_df$start_pos ] > earliest_time_alt
            
            # check duration is longer than minimum weardown duration
            meets_duration_criteria <- ( indiv_run$time[ indiv_run_rle_df$end_pos ] - indiv_run$time[ indiv_run_rle_df$start_pos ] ) > weardown_dur_alt
            
            # meets all criteria (including above-threshold)?
            indiv_run_rle_df$meets_all_criteria <- ( meets_timing_criteria * meets_duration_criteria * indiv_run_rle_df$value ) == 1
            
            if( sum( indiv_run_rle_df$meets_all_criteria ) > 0 ) {
                active_run_pos_alt <- c( active_run_pos_alt, (unique(sims_accept_alt2$run_num))[iter_alt2] )
                active_timepoint_alt <- rbind( active_timepoint_alt, c( indiv_run[ min( indiv_run_rle_df$start_pos[ indiv_run_rle_df$meets_all_criteria ] ), c( "run_num", "time", "log_mtb_act_gra_t", "t_cell_perc" ) ] ) )
            }
    }
}

active_timepoint_alt$time_rescaled <- round( active_timepoint_alt$time * ( 1000 / max(sims_accept_alt2$time) ), 0 )
active_timepoint_alt$color <- rev( colorRampPalette( brewer.pal(9,"Reds") )( 1200 )[ 21:1020 ] )[ active_timepoint_alt$time_rescaled ]

### alt #############################################

p1_alt_posts2 <- ggplot() +
    theme_bw() +
    theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
    geom_line( data=sims_accept_alt2, aes( x=time, y=log_mtb_act_gra_t, group=run_num ), colour="gray20", alpha=0.5 ) +
    geom_line( data=subset( sims_accept_alt2, run_num == accept_pos_alt[ which.min( combo_score_alt[ accept_pos_alt ] ) ] ), aes( x=time, y=log_mtb_act_gra_t ), colour="blue", size=1.5, alpha=0.5 ) +
    geom_line( data=sims_accept_alt2[ sims_accept_alt2$run_num %in% active_run_pos_alt, ], aes( x=time, y=log_mtb_act_gra_t, color=factor(run_num) ), size=1.5, alpha=0.5 ) +
    geom_point( data=active_timepoint_alt, aes( x=time, y=log_mtb_act_gra_t ), pch=13, size=3 ) +
    # geom_point( data=active_timepoint, aes( x=time, y=log_mtb_act_gra_t, colour=factor(run_num) ), pch=13, size=3 ) +
    xlab( "Time (days)" ) +
    ylab( "log10 CFU" ) +
    scale_color_manual( values=active_timepoint_alt$color ) +
    theme( legend.position="none" ) +
    ylim( c( 0, log10( mtb_max_gra + 1 ) ) )

p1_alt_posts3 <- ggplot() +
    theme_bw() +
    theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
    geom_line( data=sims_accept_alt2, aes( x=time, y=t_cell_perc, group=run_num ), colour="gray20", alpha=0.5 ) +
    geom_line( data=subset( sims_accept_alt2, run_num == accept_pos_alt[ which.min( combo_score_alt[ accept_pos_alt ] ) ] ), aes( x=time, y=t_cell_perc ), colour="blue", size=1.5, alpha=0.5 ) +
    geom_line( data=sims_accept_alt2[ sims_accept_alt2$run_num %in% active_run_pos_alt, ], aes( x=time, y=t_cell_perc, colour=factor(run_num) ), size=1.5, alpha=0.5 ) +
    geom_point( data=active_timepoint_alt, aes( x=time, y=t_cell_perc ), pch=13, size=3 ) +
    # geom_point( data=active_timepoint, aes( x=time, y=log_mtb_act_gra_t, colour=factor(run_num) ), pch=13, size=3 ) +
    xlab( "Time (days)" ) +
    ylab( "Treg %age of CD4" ) +
    ylim( c( 0, 100 ) ) +
    scale_color_manual( values=active_timepoint_alt$color ) +
    theme( legend.position="none" )


incid_alt_rle <- rle( active_timepoint_alt[ order( active_timepoint_alt$time ), ]$time )
incid_alt_rle_df <- data.frame( time=incid_alt_rle$values, incid=(cumsum(incid_alt_rle$lengths))/length(posts_alt_list) )

ann_times <- seq( from=start_time, to=end_time, by=365 )[-1]
ann_incid <- unlist( lapply( ann_times, function(x){ incid_alt_rle_df$incid[ max( which( incid_alt_rle_df$time <= x ) ) ] } ) )
ann_incid <- c( ann_incid[1], diff( ann_incid ) )       # http://stackoverflow.com/questions/21418287/functional-way-to-reverse-cumulative-sum
ann_incid_df <- data.frame( time=ann_times, incid=ann_incid )
cumul_incid_alt_rle_df <- data.frame( time=c( 0, unlist( lapply( incid_alt_rle_df$time, function(x){ rep(x,2) } ) ), end_time ), incid=c( 0, 0, unlist( lapply( rev(rev(incid_alt_rle_df$incid)), function(x){ rep(x,2) } ) ) ) )

kamat_data <- read.csv( "kamat_data.csv" )
kamat_data$cumul <- cumsum( kamat_data$Total )
cumul_kamat_data <- data.frame(
    time=c( 0, unlist( lapply( kamat_data$Time[-1], function(x){ rep(x,2) } ) )-0.5, tail( kamat_data$Time,1 )+0.5 ),
    incid=c( unlist( lapply( rev(rev(kamat_data$cumul)), function(x){ rep(x,2) } ) ) ) )

incid_alt_plot <- ggplot( ) +
    theme_bw() +
    theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
    # geom_point( data=incid_alt_rle_df, aes( x=time, y=incid ), pch=13, size=4 ) +
    geom_line( data=cumul_incid_alt_rle_df, aes( x=time/365, y=incid*100, color="Model" ), size=1.5, alpha=0.8 ) +
    geom_point( data=ann_incid_df, aes( x=time/365, y=incid*100, color="Model"), pch=15, size=3.5, alpha=0.6 ) +
    geom_point( data=ann_incid_df, aes( x=time/365, y=incid*100 ), pch=0, size=3.5, color="black", alpha=0.6 ) +
    xlab( "Time (years)" ) +
    ylab( "Incidence (% of population)" ) +
    # xlim( c( 0, end_time/365 ) ) +
    # geom_line( data=kamat_data, aes( x=Time, y=cumul ), color="red" ) +
    geom_line( data=cumul_kamat_data, aes( x=time, y=incid*100, color="Kamat" ), size=1.5, alpha=0.8 ) +
    geom_point( data=kamat_data, aes( x=Time, y=Total*100, color="Kamat" ), pch=19, size=4, alpha=0.6 ) +
    geom_point( data=kamat_data, aes( x=Time, y=Total*100 ), pch=1, size=4, color="black", alpha=0.6 ) +
    scale_x_continuous( breaks = c( 0, 2, 4, 6, 8, 10 ) ) +
    # scale_y_continuous( breaks = c( 0, 2, 4, 6, 8, 10, 12 ) ) +
    scale_y_continuous( breaks = seq( 0, round( 100*max( cumul_incid_alt_rle_df$incid, cumul_kamat_data$y ) ), 2 ) ) +
    scale_color_manual( name="Source", values=c( Model="blue", Kamat="red" ) ) +
    labs( caption = paste( "threshold=", threshold_gran_alt, ", earliest time=", earliest_time_alt, sep="" ) )

pdf( height=5, width=8, "output_fitting_alt_hists2.pdf" )
print( p1_alt_posts2 )
print( p1_alt_posts3 )
print( incid_alt_plot )
dev.off()

write.csv( sims_accept_alt2, "sims_accept_alt.csv", row.names=F )
write.csv( active_timepoint_alt, "active_timepoint_alt.csv", row.names=F )

#########################################################


