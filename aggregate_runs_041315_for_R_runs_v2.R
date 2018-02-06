####### Cluster all the runs


# file_dir <- "C:/Users/stewchang/Documents/EMOD"
# file_dir <- "Q:/TB/Intrahost_model_prototypes/reactivation_timer/philip_adaptation_032315c"
# file_dir <- "Q:/TB/Intrahost_model_prototypes/reactivation_timer/philip_adaptation_041315"
# file_dir <- "C:/Users/stewchang/Documents/EMOD/intrahost_candidate_models_111315d_versions/reduce_incidence2c5"
# file_dir <- "C:/Users/stewchang/Documents/EMOD/intrahost_candidate_models_111315d_versions/reduce_incidence2c5 - copy 170513 - seed 3 - fewer time windows - smallest priors - more runs"
# file_dir <- "C:/Users/stewchang/Documents/EMOD/intrahost_candidate_models_111315d_versions/reduce_incidence2c5 - copy 170512 - seed 1 - fewer time windows - smaller priors - more runs redo"

# setwd( "C:/Users/stewchang/Documents/EMOD/intrahost_candidate_models_111315d_versions/reduce_incidence2c5 - copy 170524 - new treg" )
setwd( "/home/ubuntu/treg07" )

#### Choose one of the following:

# file_list <- c( "time_courses.csv" )
# file_stub <- ""

file_list <- c( "time_courses_alt.csv" )
file_stub <- "_alt"

# file_list <- c( "time_courses_z.csv" )
# file_stub <- "_z"

######################

library( reshape )
# library( dtw )
library( ggplot2 )
library( RColorBrewer )

# num_clust <- 3

flynn_data <- read.csv( "Lin_Flynn_2014_main_Fig4a.csv" )

####################################

# file_list <- dir()[ grep( "test_gran_\\d", dir(), perl=T ) ]
# file_list <- rep( dir()[ grep( "test_gran", dir() ) ][1], 3 )
# file_list <- dir()[ grep( "200reps_test_gran", dir() ) ][1:2]
# file_list <- dir()[ grep( "200reps_test_gran", dir() ) ]
# file_list <- dir()[ grep( "_gran", dir() ) ]
# file_list <- dir()[ grep( "_gran", dir() ) ]

# data_raw <- read.csv( "test_gran.txt", header=T )

### http://stackoverflow.com/questions/21891841/importing-only-every-nth-row-from-a-csv-file-in-r
### Requires installation of gawk: http://sourceforge.net/projects/gnuwin32/?source=typ_redirect and installation into PATH

# write.csv( 1:1000, file='test.csv' )
# file.pipe <- pipe( "awk 'BEGIN{i=0}{i++;if (i%4==0) print $1}' < test.csv " )
# res <- read.csv( file.pipe )

big_data <- data.frame()

# # iter <- 1
# for ( iter in 1:length( file_list ) ) {
    
    # print( paste( iter, file_list[ iter ] ) )
    
    # # tmp <- readLines( pipe( 'paste( file_dir, file_list[iter], sep="/" )' ) )
    
    # # file.pipe <- pipe( paste( "gawk 'BEGIN{i=0}{if (i%500==1) print $1;i++}' < ", paste( file_dir, file_list[iter], sep="/" ) ) )
    # # file.pipe <- pipe( paste( "gawk 'BEGIN { FS = \",\" }{ if ( $1 ~ /[0-1]0$/ ) print $1, $2 }' < ", paste( file_dir, file_list[iter], sep="/" ) ) )
    # ### http://web.mit.edu/gnu/doc/html/gawk_8.html#SEC49
    # ### grab 1st and 3rd columns -- identifier and active Mtb number
    # file.pipe <- pipe( paste( "gawk 'BEGIN { FS = \",\" ; regexp = \"_0$|00$|50$\" }{ if ( $1 ~ regexp ) print $1, $3 }' < ", paste( file_dir, file_list[iter], sep="/" ) ) )      # doesn't work on work-computer -- works with gawk 4+
    
    
    # # file.pipe <- pipe( paste( "gawk 'BEGIN { FS = \",\" ; regexp = \"_0$|25$|50$|75$|00$|\" }{ if ( $1 ~ regexp ) print $1, $2 }' < ", paste( file_dir, file_list[iter], sep="/" ) ) )       # doesn't work on work-computer -- works with gawk 4+
    # # file.pipe <- pipe( paste( "gawk \"{ print $1 }\" < ", file_list[iter] ) )     # works but stripped down
    # # file.pipe <- pipe( paste( "gawk \"BEGIN { foo=1 }{ print $1 }\" < ", file_list[iter] ) )     # works but stripped down
    # # file.pipe <- pipe( paste( "gawk \"BEGIN { foo=1 }{ print $1 }\" < ", "test_gran.txt" ) )     # works but stripped down
    # # file.pipe <- pipe( paste( "gawk -F, \"{ print $1 }\" < ", "test_gran.txt" ) )     # works but stripped down
    # # file.pipe <- pipe( paste( "gawk -F, \"BEGIN{}{ print $1 }\" < ", "test_gran.txt" ) )     # works but stripped down
    # # file.pipe <- pipe( paste( "gawk -F, \"BEGIN{ regexp = 1 }{ print $1 }\" < ", "test_gran.txt" ) )     # works but stripped down
    # # file.pipe <- pipe( paste( "gawk -F, \"BEGIN{ regexp = \"50\" }{ print $1 }\" < ", "test_gran.txt" ) )     # works but stripped down
    # # file.pipe <- pipe( paste( "gawk -F, \"BEGIN{ regexp = \"50\" }{ if ( $1 ~ regexp ) print $1 }\" < ", "test_gran.txt" ) )     # works but stripped down
    
    # ###
    
    # if ( iter == 1 ) {
        # big_data <- read.csv( file.pipe, header=F, sep=" " )
        # # big_data <- read.csv( file.pipe, header=T, sep=" " )
        # # print( head( big_data ) )
        # print( tail( big_data ) )
        
    # } else {
        # # prev_max <- max( as.numeric( unlist( lapply( strsplit( levels( big_data[,1] ), split="_" ), function(x){ x[1] } ) ) ) )
        # prev_max <- max( as.numeric( unlist( lapply( strsplit( as.character( big_data[,1] ), "_" ), function(x) { x[1] } ) ) ) )
        # print( prev_max )
        
        # temp_data <- read.csv( file=file.pipe, header=F, sep=" " )
        # # print( tail( temp_data ) )
        # # # # lapply( strsplit( levels( temp_1[,1] ), split="_" ), function(x){ paste( ( as.numeric( x[1] ) + 1 ), x[2], sep="_" ) } )
        # temp_data[,1] <- unlist( lapply( strsplit( as.character( temp_data[,1] ), split="_" ), function(x){ paste( as.numeric( x[1] ) + as.numeric( prev_max ) + 1, x[2], sep="_" ) } ) )
        # print( tail( temp_data ) )
        
        # big_data <- rbind( big_data, temp_data )
        # # # print( head( big_data ) )
    # }
    
# }

if ( length( file_list ) == 1 ) {       # customize for single run from R
    big_data <- read.csv( file_list, header=T, sep="," )
    bad_runs <- c()
    for ( i in 1:length(unique(big_data$run_num)) ) {        # exclude any NA's
        if( sum( is.na( big_data[ big_data$run_num == unique(big_data$run_num)[i], "log_mtb_act_gra_t" ] ) ) > 0 ) {
            bad_runs <- c( bad_runs, unique(big_data$run_num)[i] )
        }
    }
    if ( length( bad_runs ) > 0 ) {
        big_data <- big_data[ -which( big_data$run_num %in% bad_runs ), ]
    }
    big_data$time <- round( big_data$time, 0 )
    big_data <- big_data[ which( big_data$time %in% unique(big_data$time)[ seq( 1, length(unique(big_data$time)), 1 ) ] ), ]        # Take everyt Nth time point
    big_data$id <- with( big_data, paste0( run_num, "_", time ) )
    big_data <- big_data[ , c( "id", "mtb_act_gra_t" ) ]
}


# temp_data <- read.csv( "test_gran.txt", header=F, skip=1 )
# temp_data_0 <- read.csv( "test_gran_0.txt", header=F, skip=1 )
# temp_data_1 <- read.csv( "test_gran_1.txt", header=F, skip=1 )

# temp_data[,1] <- unlist( lapply( strsplit( levels( temp_data[,1] ), split="_" ), function(x){ paste( as.numeric( x[1] ) + as.numeric( prev_max ), x[2], sep="_" ) } ) )

# data_raw <- read.csv( "test_gran.txt", header=T )
data_raw <- big_data
# colnames( data_raw ) <- c( "id", "act", "air", "qui", "spu" )
colnames( data_raw ) <- c( "id", "act" )

temp_raw <- data.frame( id=data_raw[,1] )
# temp_raw <- cbind( temp_raw, apply( data_raw[,-1], 2, function(x){ log10( x + 1 ) } ) )
temp_raw <- cbind( temp_raw, act=unlist( lapply( data_raw[,-1], function(x){ log10( x + 1 ) } ) ) )
data_raw <- temp_raw

# data_split <- within( data_raw, X <- data.frame( do.call( 'rbind', strsplit( as.character( X ), "_", fixed=TRUE ) ) ) )
X <- data.frame( do.call( 'rbind', strsplit( as.character( data_raw$id ), "_", fixed=TRUE ) ) )      # http://stackoverflow.com/questions/7069076/split-column-at-delimiter-in-data-frame
data_processed <- data.frame( indivID=X$X1, tpoint=X$X2, data_raw )
data_processed <- data_processed[ , -grep( "id", colnames( data_processed ) ) ]

data_melted <- melt( data_processed, id=c( "indivID", "tpoint" ) )
#data_melted$tpoint <- as.numeric( levels( data_melted$tpoint ) )
#temp <- data_melted[ intersect( which( data_melted$indivID == "11" ), which( data_melted$variable == "act" ) ), ]
#plot( temp$value )

#data_melted[ which( data_melted$variable == "act" ), ]

### make time-series ready
data_processed2 = data.frame( matrix( nrow = 1, ncol = length( unique( data_melted$tpoint ) ) ) )
row_names <- list()
indivIDs <- as.character( unique( data_melted$indivID ) )
for ( i in 1:length( indivIDs ) ) {
    row_names <- append( row_names, indivIDs[i] )
    temp <- data_melted$value[ which( data_melted$indivID == indivIDs[i] ) ]
    # print( temp )
    data_processed2 <- rbind( data_processed2, temp )
}
data_processed2 <- data_processed2[ -1, ]
row_names <- unlist( row_names )
rownames( data_processed2 ) <- row_names
colnames( data_processed2 ) <- as.character( data_melted$tpoint[ which( data_melted$tpoint == 0 )[1]:(which( data_melted$tpoint == 0 )[2]-1) ] )
# data_processed3 <- 10**data_processed2-1

#############

distMatrix2 <- dist(data_processed2, method="euclidean")
# distMatrix3 <- dist(data_processed3, method="euclidean")

num_clust <- 3

# pdf( "clustering_output_euclidean_log10.pdf", height=6, width=9 )
pdf( paste( "clustering_output_euclidean_log10", file_stub, ".pdf", sep="" ), height=6, width=9 )

hc2 <- hclust( distMatrix2, method="average" )
plot( hc2, main="Clustering by Euclidean distance" )                                             # labels are row numbers, can be used to check
rect.hclust( hc2, k=num_clust, border="red" )
plot( hc2, main="Clustering by Euclidean distance", labels=F )                                             # labels are row numbers, can be used to check
hclust2_out <- rect.hclust( hc2, k=num_clust, border="red" )
#myhc <- cutree( hc, h=100 )                                   # hc plot shows where to apply the cut
myhc2 <- cutree( hc2, k=num_clust )                                   # hc plot shows where to apply the cut

raw_widths <- unlist( append( 0, ( lapply( hclust2_out, function(x){ length(x) } ) ) ) )
xpos <- cumsum( raw_widths )[ -length( cumsum( raw_widths ) ) ] + diff( cumsum( raw_widths ) ) / 2
clust_labels <- unlist( lapply( hclust2_out, function(x){ names( table(myhc2)[ match( length(x), table(myhc2) ) ] ) } ) )
# text( x=xpos, y=-1.75*sd(hc2$height), labels=clust_labels, col="red" )
text( x=xpos, y=-max(hc2$height)/8, labels=clust_labels, col="red" )

for ( j in 1:length( unique( myhc2 ) ) ) {
    print( j )
    
    group <- which( myhc2 == unique( myhc2 )[ j ] )
    # print( group )
    
    group1_val <- data_processed2[ group[1], ]
    print( group1_val )
    plot( x=as.numeric( colnames( data_processed2 ) ), y=as.numeric( group1_val ), type="l", ylim=c(0,5), xlab="Time (d)", ylab="Log10( Number of active Mtb + 1)", main=paste( "Cluster ", unique(myhc2)[j], " (n=", length(group), ")", sep="" ) )
    
    for ( i in names( group )[-1] ) {
        points( x=as.numeric( colnames( data_processed2 ) ), y=data_processed2[i,], type="l" )
    }
}

data_processed2_melted <- melt( t( data_processed2 ) )
colnames( data_processed2_melted ) <- c( "timepoint", "id", "value" )
cols <- brewer.pal( num_clust, "Set2" )

for ( j in unique( myhc2 ) ) {
    print( j )
    group1 <- names( which( myhc2 == j ) )
    group1_data <- data_processed2_melted[ which( as.character( data_processed2_melted$id ) %in% group1 ), ]
    
    # group1_data_byTime <- data.frame( timepoint2 = group1_data$time[ group1_data$id == group1_data$id[1] ], means = colMeans( cast( group1_data, id ~ timepoint ) ) )
    group1_data_recast <- cast( group1_data, id ~ timepoint )
    group1_data_byTime <- data.frame( timepoint2 = group1_data$time[ group1_data$id == group1_data$id[1] ],
                                      means = apply( cast( group1_data, id ~ timepoint ), 2, function(x){ quantile( x, probs = c( 0.5 ) ) } ),      # note, actually now using median
                                      lower = apply( cast( group1_data, id ~ timepoint ), 2, function(x){ quantile( x, probs = c( 0.05 ) ) } ),
                                      upper = apply( cast( group1_data, id ~ timepoint ), 2, function(x){ quantile( x, probs = c( 0.95 ) ) } ) )
                                      # means = colMeans( group1_data_recast ),
                                      # lower = apply( cast( group1_data, id ~ timepoint ), 2, function(x){ mean(x)-( qnorm(0.975)*sd(x)/sqrt(length(x)) ) } ),
                                      # upper = apply( cast( group1_data, id ~ timepoint ), 2, function(x){ mean(x)+( qnorm(0.975)*sd(x)/sqrt(length(x)) ) } ) )    
    
    # p1 <- ggplot( data=group1_data, aes( x=timepoint, y=value, group=id ) ) +
    p1 <- ggplot( ) +
        # geom_line( aes(alpha=0.05), color=cols[grep(j,unique(myhc2))] ) +
        geom_line( data=group1_data, aes( x=timepoint, y=value, group=id, alpha=0.05 ), color=cols[grep(j,unique(myhc2))] ) +
        theme_bw() +
        theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
        geom_point( data=data.frame( Time=flynn_data$Time, log_CFU=log10( flynn_data$CFU ) ), aes( x=Time, y=log_CFU ) ) +
        # geom_line( data=group1_data_byTime, aes( x=timepoint2, y=means ), color="gray20" ) +
        # geom_line( data=group1_data_byTime, aes( x=timepoint2, y=lower ), color="gray20" ) +
        # geom_line( data=group1_data_byTime, aes( x=timepoint2, y=upper ), color="gray20" ) +
        ylim( c(0,5) ) +
        xlab( "Time (days)" ) +
        ylab( "log10 CFU" ) +
        ggtitle( paste( "Cluster ", unique(myhc2)[j], " (n=", length(group1), ")", sep="" ) ) +
        theme( legend.position="none" )
    print( p1 )
}

# dev.off()

data_processed2_early <- data_processed2[ , which( as.numeric( colnames( data_processed2 ) ) <= 500 ) ]
for ( j in 1:length( levels( factor( myhc2 ) ) ) ) {
    group <- which( myhc2 == levels( factor( myhc2 ) )[ j ] )
    group1_val <- data_processed2_early[ group[1], ]
    plot( x=as.numeric( colnames( data_processed2_early ) ), y=as.numeric( group1_val ), type="l", ylim=c(0,5), xlab="Time (d)", ylab="Log10( Number of active Mtb + 1)", main=paste( "Cluster ", unique(myhc2)[j], " (n=", length(group), ")", sep="" ) )
    
    for ( i in names( group )[-1] ) {
        points( x=as.numeric( colnames( data_processed2_early ) ), y=data_processed2_early[i,], type="l" )
    }
}

for ( j in unique( myhc2 ) ) {
    print( j )
    group1 <- names( which( myhc2 == j ) )
    group1_data <- data_processed2_melted[ intersect( which( as.character( data_processed2_melted$id ) %in% group1 ), which( data_processed2_melted$timepoint <= 500 ) ), ]         # times less than 500 d
    
    group1_data_byTime <- data.frame( timepoint2 = group1_data$time[ group1_data$id == group1_data$id[1] ],
                                      means = apply( cast( group1_data, id ~ timepoint ), 2, function(x){ quantile( x, probs = c( 0.5 ) ) } ),      # note, actually now using median
                                      lower = apply( cast( group1_data, id ~ timepoint ), 2, function(x){ quantile( x, probs = c( 0.05 ) ) } ),
                                      upper = apply( cast( group1_data, id ~ timepoint ), 2, function(x){ quantile( x, probs = c( 0.95 ) ) } ) )
                                      # means = colMeans( cast( group1_data, id ~ timepoint ) ) )
    
    # p1 <- ggplot( group1_data, aes( x=timepoint, y=value, group=id ) ) +
    p1 <- ggplot( ) +
        # geom_line( aes(alpha=0.05), color=cols[grep(j,unique(myhc2))] ) +
        geom_line( data=group1_data, aes( x=timepoint, y=value, group=id, alpha=0.05 ), color=cols[grep(j,unique(myhc2))] ) +
        theme_bw() +
        theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
        geom_point( data=data.frame( Time=flynn_data$Time, log_CFU=log10( flynn_data$CFU ) ), aes( x=Time, y=log_CFU ) ) +
        # geom_line( data=group1_data_byTime, aes( x=timepoint2, y=means ), color="gray20" ) +
        # geom_line( data=group1_data_byTime, aes( x=timepoint2, y=lower ), color="gray20" ) +
        # geom_line( data=group1_data_byTime, aes( x=timepoint2, y=upper ), color="gray20" ) +
        ylim( c(0,5) ) +
        xlab( "Time (days)" ) +
        ylab( "log10 CFU" ) +
        ggtitle( paste( "Cluster ", unique(myhc2)[j], " (n=", length(group1), ")", sep="" ) ) +
        theme( legend.position="none" )
    print( p1 )
}

dev.off()

# write.csv( myhc2, file="clustering_output_euclidean_log10_clusters.csv" )
write.csv( myhc2, file=paste( "clustering_output_euclidean_log10_clusters", file_stub, ".csv", sep="" ) )


############

# distMatrix1 <- dist(data_processed2, method="DTW")

# num_clust <- 3

# # pdf( "clustering_output_dtw_log10.pdf", height=6, width=9 )
# pdf( paste( "clustering_output_dtw_log10", file_stub, ".pdf", sep="" ), height=6, width=9 )

# hc1 <- hclust( distMatrix1, method="average" )
# plot( hc1, main="Clustering by DTW" )                                             # labels are row numbers, can be used to check
# rect.hclust( hc1, k=num_clust, border="red" )
# plot( hc1, main="Clustering by DTW", labels=F )                                             # labels are row numbers, can be used to check
# hclust1_out <- rect.hclust( hc1, k=num_clust, border="red" )
# #myhc <- cutree( hc, h=100 )                                   # hc plot shows where to apply the cut
# myhc1 <- cutree( hc1, k=num_clust )                                   # hc plot shows where to apply the cut

# raw_widths <- unlist( append( 0, ( lapply( hclust1_out, function(x){ length(x) } ) ) ) )
# xpos <- cumsum( raw_widths )[ -length( cumsum( raw_widths ) ) ] + diff( cumsum( raw_widths ) ) / 2
# clust_labels <- unlist( lapply( hclust1_out, function(x){ names( table(myhc1)[ match( length(x), table(myhc1) ) ] ) } ) )
# # text( x=xpos, y=-3*sd(hc1$height), labels=clust_labels, col="red" )
# text( x=xpos, y=-max(hc1$height)/8, labels=clust_labels, col="red" )

# # for ( j in 1:length( levels( factor( myhc ) ) ) ) {
# for ( j in 1:length( unique( myhc1 ) ) ) {
    # print( j )
    
    # # group <- which( myhc == levels( factor( myhc ) )[ j ] )
    # group <- which( myhc1 == unique( myhc1 )[ j ] )
    # print( group )
    
    # # group1_val <- data_processed3_small[ group[1], ]
    # group1_val <- data_processed2[ group[1], ]
    # print( group1_val )
    # #ggplot( data=subset( data_melted, indivID %in% group_indivID & variable == "act" ) ) +
    # #    geom_point( aes( x=tpoint, y=value ) )
    # # plot( x=as.numeric( colnames( data_processed3_small ) ), y=as.numeric( group1_val ), type="l", ylim=c(0,5), xlab="Time (d)", ylab="Log10( Number of active Mtb + 1)", main=paste("Cluster ",j) )
    # plot( x=as.numeric( colnames( data_processed2 ) ), y=as.numeric( group1_val ), type="l", ylim=c(0,5), xlab="Time (d)", ylab="Log10( Number of active Mtb + 1)", main=paste( "Cluster ", unique(myhc1)[j], " (n=", length(group), ")", sep="" ) )
    
    # for ( i in names( group )[-1] ) {
        # # points( x=as.numeric( colnames( data_processed3_small ) ), y=data_processed3_small[i,], type="l" )
        # points( x=as.numeric( colnames( data_processed2 ) ), y=data_processed2[i,], type="l" )
    # }
# }

# data_processed2_melted <- melt( t( data_processed2 ) )
# colnames( data_processed2_melted ) <- c( "timepoint", "id", "value" )
# cols <- brewer.pal( num_clust, "Set2" )

# for ( j in unique( myhc1 ) ) {
    # print( j )
    # group1 <- names( which( myhc1 == j ) )
    # group1_data <- data_processed2_melted[ which( as.character( data_processed2_melted$id ) %in% group1 ), ]
    
    # group1_data_recast <- cast( group1_data, id ~ timepoint )
    # group1_data_byTime <- data.frame( timepoint2 = group1_data$time[ group1_data$id == group1_data$id[1] ],
                                      # means = apply( cast( group1_data, id ~ timepoint ), 2, function(x){ quantile( x, probs = c( 0.5 ) ) } ),      # note, actually now using median
                                      # lower = apply( cast( group1_data, id ~ timepoint ), 2, function(x){ quantile( x, probs = c( 0.05 ) ) } ),
                                      # upper = apply( cast( group1_data, id ~ timepoint ), 2, function(x){ quantile( x, probs = c( 0.95 ) ) } ) )
                                      # # means = colMeans( group1_data_recast ),
                                      # # lower = apply( cast( group1_data, id ~ timepoint ), 2, function(x){ mean(x)-( qnorm(0.975)*sd(x)/sqrt(length(x)) ) } ),
                                      # # upper = apply( cast( group1_data, id ~ timepoint ), 2, function(x){ mean(x)+( qnorm(0.975)*sd(x)/sqrt(length(x)) ) } ) )
    
    # # p1 <- ggplot( data=group1_data, aes( x=timepoint, y=value, group=id ) ) +
    # p1 <- ggplot( ) +
        # # geom_line( aes(alpha=0.05), color=cols[grep(j,unique(myhc2))] ) +
        # geom_line( data=group1_data, aes( x=timepoint, y=value, group=id, alpha=0.05 ), color=cols[grep(j,unique(myhc1))] ) +
        # theme_bw() +
        # theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
        # geom_point( data=data.frame( Time=flynn_data$Time, log_CFU=log10( flynn_data$CFU ) ), aes( x=Time, y=log_CFU ) ) +
        # # geom_line( data=group1_data_byTime, aes( x=timepoint2, y=means ), color="gray20" ) +
        # # geom_line( data=group1_data_byTime, aes( x=timepoint2, y=lower ), color="gray20" ) +
        # # geom_line( data=group1_data_byTime, aes( x=timepoint2, y=upper ), color="gray20" ) +
        # ylim( c(0,5) ) +
        # xlab( "Time (days)" ) +
        # ylab( "log10 CFU" ) +
        # ggtitle( paste( "Cluster ", unique(myhc1)[j], " (n=", length(group1), ")", sep="" ) ) +
        # theme( legend.position="none" )
    # print( p1 )
# }

# # dev.off()

# data_processed2_early <- data_processed2[ , which( as.numeric( colnames( data_processed2 ) ) <= 500 ) ]
# for ( j in 1:length( levels( factor( myhc1 ) ) ) ) {
    # group <- which( myhc1 == levels( factor( myhc1 ) )[ j ] )
    # group1_val <- data_processed2_early[ group[1], ]
    # plot( x=as.numeric( colnames( data_processed2_early ) ), y=as.numeric( group1_val ), type="l", ylim=c(0,5), xlab="Time (d)", ylab="Log10( Number of active Mtb + 1)", main=paste( "Cluster ", unique(myhc1)[j], " (n=", length(group), ")", sep="" ) )
    
    # for ( i in names( group )[-1] ) {
        # points( x=as.numeric( colnames( data_processed2_early ) ), y=data_processed2_early[i,], type="l" )
    # }
# }

# for ( j in unique( myhc1 ) ) {
    # print( j )
    # group1 <- names( which( myhc1 == j ) )
    # group1_data <- data_processed2_melted[ intersect( which( as.character( data_processed2_melted$id ) %in% group1 ), which( data_processed2_melted$timepoint <= 500 ) ), ]         # times less than 500 d
    
    # group1_data_byTime <- data.frame( timepoint2 = group1_data$time[ group1_data$id == group1_data$id[1] ],
                                      # means = apply( cast( group1_data, id ~ timepoint ), 2, function(x){ quantile( x, probs = c( 0.5 ) ) } ),      # note, actually now using median
                                      # lower = apply( cast( group1_data, id ~ timepoint ), 2, function(x){ quantile( x, probs = c( 0.05 ) ) } ),
                                      # upper = apply( cast( group1_data, id ~ timepoint ), 2, function(x){ quantile( x, probs = c( 0.95 ) ) } ) )
                                      # # means = colMeans( cast( group1_data, id ~ timepoint ) ) )
    
    # # p1 <- ggplot( group1_data, aes( x=timepoint, y=value, group=id ) ) +
    # p1 <- ggplot( ) +
        # # geom_line( aes(alpha=0.05), color=cols[grep(j,unique(myhc2))] ) +
        # geom_line( data=group1_data, aes( x=timepoint, y=value, group=id, alpha=0.05 ), color=cols[grep(j,unique(myhc1))] ) +
        # theme_bw() +
        # theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
        # geom_point( data=data.frame( Time=flynn_data$Time, log_CFU=log10( flynn_data$CFU ) ), aes( x=Time, y=log_CFU ) ) +
        # # geom_line( data=group1_data_byTime, aes( x=timepoint2, y=means ), color="gray20" ) +
        # # geom_line( data=group1_data_byTime, aes( x=timepoint2, y=lower ), color="gray20" ) +
        # # geom_line( data=group1_data_byTime, aes( x=timepoint2, y=upper ), color="gray20" ) +
        # ylim( c(0,5) ) +
        # xlab( "Time (days)" ) +
        # ylab( "log10 CFU" ) +
        # ggtitle( paste( "Cluster ", unique(myhc1)[j], " (n=", length(group1), ")", sep="" ) ) +
        # theme( legend.position="none" )
    # print( p1 )
# }

# dev.off()

# # write.csv( myhc1, file="clustering_output_dtw_log10_clusters.csv" )
# write.csv( myhc1, file=paste( "clustering_output_dtw_log10_clusters", file_stub, ".csv", sep="" ) )

