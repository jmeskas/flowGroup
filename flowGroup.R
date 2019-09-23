# Written by Justin Meskas
# Last updated on September 2019
# alpha version 0.99.5

flowGroup <- function(files, xMark, yMark, experiment_name = NULL, plot_groups = F, HOG_cex=0.2,
                      grid_width = 200, partitions_flowType = 4, vuc = 1, vuc_d = 1, hog = 1, flwT = 1,
                      vert_line = NULL, horiz_line = NULL, xlim=NA, ylim=NA, load_files = F, Labels=NULL,
                      DoOverlay = F, DoContour = F, boundaries=NULL, k2run=NULL, verbose = F ){

    start <- Sys.time()

    switch(class(files),
           flowFrame = {flowFrameOrPathnames <- "flowFrame"},
           flowSet   = {flowFrameOrPathnames <- "flowFrame"},
           list      = {ifelse(class(files[[1]]) == "flowFrame",
                               flowFrameOrPathnames <- "flowFrame",
                               flowFrameOrPathnames <- "Pathnames")},
           character = {flowFrameOrPathnames <- "Pathnames"})

    if (load_files == T && flowFrameOrPathnames == "Pathnames"){
        files <- unlist(files)

        files <- lapply(files, function(x){read.FCS(filename = paste0(x))})
    }

    date_time <- strftime(Sys.time(), "%y%m%d_%H%M")

    no_cores <- detectCores() - 8
    if (no_cores <= 0){
        no_cores <- 1
    }
    registerDoMC(no_cores)

    if ( is.null(experiment_name)){ # need a folder name to store all the figures
        experiment_name <- date_time
        if(verbose==T){ cat("Experiment name set to: ", date_time,"\n",sep="") }
    }
    suppressWarnings( dir.create(experiment_name, recursive = T) )
    suppressWarnings( dir.create(paste0(experiment_name,"/images/"), recursive = T) )
    suppressWarnings( dir.create(paste0(experiment_name,"/Derivatives/"), recursive = T) )

    num_of_files <- length(files)
    
    ########################################################################
    # VUC
    # source("code/flowGroup_functions.R")
    if (vuc_d != 0 || vuc != 0) {
        ########################################################################
        # density representation

        start_kde2d <- Sys.time()
        kde2d_res <- foreach(q1 = 1:num_of_files ) %dopar% {
            f_temp <- load_flowFrame(files=files, q1=q1, flowFrameOrPathnames = flowFrameOrPathnames)
            f1 <- kde2d(f_temp@exprs[,xMark], f_temp@exprs[,yMark], n=grid_width, lims=c(xlim,ylim))
            return(f1) # this returns roughly the same value. On one test it can a +- 1.6% difference.
        }
    }
        ####
        # calculate distant matrix
    if (vuc != 0) {
        volume_perc_under_both_curves <- foreach(q1 = 1:num_of_files, .combine = rbind, .maxcombine = num_of_files, .multicombine = T ) %dopar% {
            f1 <- kde2d_res[[q1]]
            ratio_value <- vector(length = num_of_files)
            for ( q2 in 1:num_of_files){
                if (q2 >= q1) {
                    ratio_value[q2] <- 0
                } else {
                    f2 <- kde2d_res[[q2]]
                    trap_value  <- sapply(1:grid_width^2, function(x) { abs(  f1$z[x] - f2$z[x]) })
                    max_value   <- sapply(1:grid_width^2, function(x) { max(c(f1$z[x],  f2$z[x])) })
                    ratio_value[q2] <- sum(trap_value)/sum(max_value)
                }
            }
            return(ratio_value)
        }
        row.names(volume_perc_under_both_curves) <- 1:num_of_files
        # distance is not euclidean or manhattan, it is my own custom distance.
        dist_vuc <- stats::as.dist(volume_perc_under_both_curves) # I use as.dist because I want my values to be in the distance matrix. I dont want my points wo be thought of as coordinates and for dist() to calculate something I dont want / is incorrect.

        # hclust and dendrogram
        hc_vuc <- hclust(dist_vuc)
        png(paste0(experiment_name, "/", date_time,"_VUC_Dendrogram.png"),width=1500,height=1000)
            plot(hc_vuc, main="Volume under the curve")
        dev.off()

        dist_vuc_norm <- dist_vuc/max(dist_vuc)

        if(verbose==T){ cat("Time for the volume under both curves method: ", TimeOutput(start_kde2d),"\n",sep="") }

    } else {
        dist_vuc <- stats::dist(matrix(0,num_of_files,num_of_files))
        dist_vuc_norm <- stats::dist(matrix(0,num_of_files,num_of_files))
        if(verbose==T){ cat("Skipped the volume under both curves method.\n",sep="") }
    }



    ########################################################################
    # VUC derivative
    if (vuc_d != 0) {

        start_kde2d <- Sys.time()
        derivative <- foreach(q1 = 1:num_of_files ) %dopar% {
            derivative = grad(h=kde2d_res[[q1]]$z, x=kde2d_res[[q1]]$x, y=kde2d_res[[q1]]$y)
            return(derivative)
        }
        doubleDerivative <- foreach(q1 = 1:num_of_files ) %dopar% {
            doubleDerivative <- sqrt(derivative[[q1]]$gx^2+derivative[[q1]]$gy^2)
            return(doubleDerivative)
        }
        Do_plots <- T
        if (Do_plots == T){
            figs_temp <- foreach(q2 = 1:num_of_files, .combine = list, .maxcombine = length(files), .multicombine = T) %dopar% {

                CairoPNG ( file = paste0(experiment_name, "/Derivatives/",str_pad(q2, 3, pad = "0"),".png"), width = 2000, height = 2000)
                    par(mfrow=c(2,2), mar=c(3, 3, 1, 1), mgp=c(2, 0.7, 0))
                    image2D(kde2d_res[[q2]]$z,      xlab="",ylab="",frame=FALSE,axes=F, main="")
                    image2D(derivative[[q2]]$gx,    xlab="",ylab="",frame=FALSE,axes=F, main="")
                    image2D(derivative[[q2]]$gy,    xlab="",ylab="",frame=FALSE,axes=F, main="")
                    image2D(doubleDerivative[[q2]], xlab="",ylab="",frame=FALSE,axes=F, main="")
                dev.off ( )
                return(1)
            }
        }

        ####
        # calculate distant matrix

        volume_perc_under_both_curves_dervative <- foreach(q1 = 1:num_of_files, .combine = rbind, .maxcombine = num_of_files, .multicombine = T ) %dopar% {
            f1 <- doubleDerivative[[q1]]
            ratio_value <- vector(length = num_of_files)
            for ( q2 in 1:num_of_files){
                if (q2 >= q1) {
                    ratio_value[q2] <- 0
                } else {
                    f2 <- doubleDerivative[[q2]]
                    trap_value  <- sapply(1:grid_width^2, function(x) { abs(  f1[[x]] - f2[[x]]) })
                    max_value   <- sapply(1:grid_width^2, function(x) { max(c(f1[[x]],  f2[[x]])) })
                    ratio_value[q2] <- sum(trap_value)/sum(max_value)
                }
            }
            return(ratio_value)
        }
        row.names(volume_perc_under_both_curves_dervative) <- 1:num_of_files
        # distance is not euclidean or manhattan, it is my own custom distance.
        dist_vuc_derv <- stats::as.dist(volume_perc_under_both_curves_dervative) # I use as.dist because I want my values to be in the distance matrix. I dont want my points wo be thought of as coordinates and for dist() to calculate something I dont want / is incorrect.

        # hclust and dendrogram
        hc_vuc <- hclust(dist_vuc_derv)
        png(paste0(experiment_name, "/", date_time,"_VUC_Derv_Dendrogram.png"),width=1500,height=1000)
            plot(hc_vuc, main="Volume under the curve of the derivatives")
        dev.off()

        dist_vuc_derv_norm <- dist_vuc_derv/max(dist_vuc_derv)

        if(verbose==T){ cat("Time for the volume under both curves of the derivatives method: ", TimeOutput(start_kde2d),"\n",sep="") }

    } else {
        dist_vuc_derv <- stats::dist(matrix(0,num_of_files,num_of_files))
        dist_vuc_derv_norm <- stats::dist(matrix(0,num_of_files,num_of_files))
        if(verbose==T){ cat("Skipped the volume under both curves of the derivative method.\n",sep="") }
    }

    ########################################################################
    # HOG plotting

    if (hog != 0) {
        start_HOG <- Sys.time()
        Do_plots <- T
        if (Do_plots == T){
            figs_temp <- foreach(q1 = 1:num_of_files, .combine = list, .maxcombine = length(files), .multicombine = T) %dopar% {

                CairoPNG ( file = paste0(experiment_name, "/images/",str_pad(q1, 3, pad = "0"),".png"), width = 4000, height = 4000)
                    par(mar=c(1,1,1,1)) # margins
                    f_temp <- load_flowFrame(files=files, q1=q1, flowFrameOrPathnames = flowFrameOrPathnames)
                    plot(f_temp@exprs[,xMark],f_temp@exprs[,yMark],xlab="",ylab="",frame=FALSE,axes=F,
                                 cex.main=2, cex.lab=2, cex.axis=2, cex=HOG_cex, pch=19, main="", xlim=xlim, ylim=ylim)
                dev.off ( )
            }
        }
        # sink("myfile", append=FALSE, split=FALSE) # want to hide the output of HOG_apply
        res_hog <- HOG_apply(paste0(experiment_name, "/images/"), cells =  100, orientations = 12, threads = no_cores)
        # sink()
        # if (file.exists("myfile")) file.remove("myfile")

        dist_hog <- stats::dist(res_hog$hog, method = "manhattan") # with cells being 100 and orientations being 12, there are 120000 data values for each file/picutre. Manhattan is better than euclidean for large dimensions.
        hc_hog <- hclust(dist_hog, method = "complete") # default
        png(paste0(experiment_name, "/", date_time,"_HOG_Dendrogram_100_4k.png"),width=1500,height=1000)
            plot(hc_hog, main="HOG", xlim=c(0,4), ylim=c(0,4))
        dev.off()

        dist_hog_norm <- dist_hog/max(dist_hog)

        if(verbose==T){ cat("Time for the HOG method: ", TimeOutput(start_HOG),"\n",sep="") }

    } else {
        dist_hog <- stats::dist(matrix(0,num_of_files,num_of_files))
        dist_hog_norm <- stats::dist(matrix(0,num_of_files,num_of_files))
        if(verbose==T){ cat("Skipped the HOG method.\n",sep="") }
    }


    ########################################################################
    # flowType partitioning
    if (flwT != 0) {
        start_flowType <- Sys.time()
        flowType_res <- foreach(q1 = 1:num_of_files ) %dopar% {
            f_temp <- load_flowFrame(files=files, q1=q1, flowFrameOrPathnames = flowFrameOrPathnames)
            if(!is.na(xlim) && !is.na(ylim)){
                xThres <- seq(from = xlim[1], to = xlim[2], length.out = partitions_flowType+1)
                yThres <- seq(from = ylim[1], to = ylim[2], length.out = partitions_flowType+1)
            } else {
                xThres <- seq(from = min(f_temp@exprs[,xMark]), to = max(f_temp@exprs[,xMark]), length.out = partitions_flowType+1)
                yThres <- seq(from = min(f_temp@exprs[,yMark]), to = max(f_temp@exprs[,yMark]), length.out = partitions_flowType+1)
            }
            xThres <- xThres[-c(1,length(xThres))] # do not want partitions of nothing on the borders
            yThres <- yThres[-c(1,length(yThres))] # do not want partitions of nothing on the borders
            res <- flowType(Frame = f_temp, PropMarkers = c(xMark, yMark), MFIMarkers = c(xMark, yMark), Methods = "Thresholds",
                            MarkerNames = c(xMark, yMark), Thresholds = list(xThres, yThres), MemLimit = 10, verbose = T,
                            PartitionsPerMarker = c(length(xThres)+1,length(yThres)+1) );
            rem.ind <- grep(x = res@PhenoCodes, pattern = "0")
            CellFreqs <- res@CellFreqs[-rem.ind]/res@CellFreqs[1]
            return(CellFreqs)
        }

        flowType_matrix <- t(cbind(sapply(1:length(flowType_res), function(x){as.vector(flowType_res[[x]])})))

        dist_flwT <- stats::dist(flowType_matrix, method = "manhattan") # choose manhattan over euclidean because there are 16 dimensions, and manhattan performs better with large amount of dimensions.
        hc_flwT <- hclust(dist_flwT)

        png(paste0(experiment_name, "/", date_time,"_flowtype_Dendrogram.png"),width=1500,height=1000)
            plot(hc_flwT, main="flowType partitions")
        dev.off()

        dist_flwT_norm <- dist_flwT/max(dist_flwT)

        if(verbose==T){ cat("Time for the flowType method: ", TimeOutput(start_flowType),"\n",sep="") }

    } else {
        dist_flwT <- stats::dist(matrix(0,num_of_files,num_of_files))
        dist_flwT_norm <- stats::dist(matrix(0,num_of_files,num_of_files))
        if(verbose==T){ cat("Skipped the flowType method.\n",sep="") }
    }
    ########################################################################
    # combine 4 versions

    dist_all <- vuc*dist_vuc_norm + vuc_d*dist_vuc_derv_norm + hog*dist_hog_norm + flwT*dist_flwT_norm

    # densities <- flowOutlierDensities(files, channel=7)
    # rank_matrix <- flowOutlierRankMatrix(densities)
    # dist_all <- stats::as.dist(t(rank_matrix))

    # hclust_group <- c("complete", "ward.D", "ward.D2", "centroid", "average", "median", "mcquitty", "single")
    # hclust_group <- hclust_group[1]

    hclust_group <- "complete"
    # date_time <- "181114"

    suppressWarnings( dir.create(paste0(experiment_name,"/", hclust_group), recursive = T) )

    dist_all_norm <- (dist_all-min(dist_all))/(max(dist_all)-min(dist_all))

    hc_all <- hclust(dist_all_norm, method = hclust_group)

    png(paste0(experiment_name, "/", hclust_group, "/", date_time,"_combined_Dendrogram_", hclust_group, ".png"),width=1500,height=1000)
        plot(hc_all, main="combined")
    dev.off()

    # if(verbose==T){ cat("Time for hclust: ", TimeOutput(start_hclust),"\n",sep="") }

    ########################################################################
    # best k values

    k_count <- NULL
    seq_t <- seq(min(hc_all$height), max(hc_all$height), length.out = 100)
    for(t1 in 1:100){
        k_count[t1] <- length(unique(cutree(hc_all, h=seq_t[t1])))
    }

    find_best_k <- sort(table(k_count),decreasing=TRUE)
    ind_rem <- sort(unique(c(which(as.numeric(names(find_best_k)) > 50), which(find_best_k <= 1)))) # dont want clusters of size 1 or k too large where the largeness is overweighting the score
    find_best_k <- find_best_k[-ind_rem]
    temp <- as.numeric(names(find_best_k))
    temp[temp>=15] <- 15
    scores <- temp*find_best_k # random scoring system. I want a medium sized k with many points. (i.e. large vertical gab on dendrogram)
    best_k <- as.numeric(names(which(scores >= 0.85*max(scores))))

    # make the table nice enough to display
    scores_vis <- t(as.matrix(scores))
    scores_vis <- rbind(scores_vis, scores_vis)
    scores_vis[1,] <- as.numeric(colnames(scores_vis))
    row.names(scores_vis) <- c("k","score")
    colnames(scores_vis) <- rep("", ncol(scores_vis))
    scores_vis <- scores_vis[,sort(scores_vis[2,], index.return=T, decreasing = T)$ix]

    if(verbose==T){ print(scores_vis) }
    # if(verbose==T){ print(paste0("best_k is ", best_k)) }

    ###############################################################
    # plot groups

    set.seed(4342)
    res.tsne <- Rtsne(dist_all_norm, theta = 0.0, perplexity = 1)
    prin_comp <- prcomp(dist_all_norm) # PCA
    return_cutree_order <- list()
    order_of_plotting_total <- list()

    index <- which(scores_vis[2,] >= 0.5*scores_vis[2,1])

    k_to_run <- sort(c(unique(c(scores_vis[1,index]))))

    if (is.null(k2run)){
        # k2run <- k_to_run[1]
        k2run <- best_k[1]
    }
    
    k_to_run <- k2run

    # k_to_run <- best_k
    for ( n_of_groups in k_to_run){
    # for ( n_of_groups in best_k){

        start_groups <- Sys.time()
        cutree_order <- cutree(hc_all,k=n_of_groups)
        names(cutree_order) <- NULL
        printing_order <- unique(cutree_order[hc_all$order])
        group_ind <- cutree_order[hc_all$order]
        width <- sort(table(group_ind),decreasing=TRUE)[1]

        order_of_plotting <- lapply(1:n_of_groups, function(x){temp <- hc_all$order[which(group_ind==printing_order[x])]
                                          c(temp,rep(num_of_files+1,width-length(temp)))})
        order_of_plotting_total[[length(order_of_plotting_total)+1]] <- order_of_plotting

        colours <- as.numeric(sapply(1:num_of_files, function(x) { names(find_number(x, order_of_plotting)) }))
        colours <- rainbow(n_of_groups+2)[colours]

        # width <- min(width, 30) # only plot 30 columns, avoids R crashing.

        return_cutree_order[[length(return_cutree_order)+1]] <- cutree_order

        # tSNE plots
        png(paste0(experiment_name, "/", hclust_group, "/tsne_", date_time,"_", n_of_groups, "_", hclust_group, ".png"),width=1000,height=1000,pointsize=18)
            plot(res.tsne$Y, col=1:num_of_files, pch='.', xlab='', ylab='')
            # text(res.tsne$Y, labels=c(1:num_of_files), cex=0.8, col=return_cutree_order[[which(n_of_groups == k_to_run)]])
            text(res.tsne$Y, labels=c(1:num_of_files), cex=0.8, col=colours)
        dev.off()

        # PCA plots
        png(paste0(experiment_name, "/", hclust_group, "/pca_", date_time,"_", n_of_groups, "_", hclust_group, ".png"),width=1600,height=1600,pointsize=18)
        par(mfrow=c(2,2), mar=c(3, 3, 1, 1), mgp=c(2, 0.7, 0))
            plot(prin_comp$x[,c(1,2)], pch='.', col=return_cutree_order[[which(n_of_groups == k_to_run)]], xlab='', ylab='')
            # text(prin_comp$x[,c(1,2)], labels=c(1:num_of_files), cex=0.8, col=return_cutree_order[[which(n_of_groups == k_to_run)]])
            text(prin_comp$x[,c(1,2)], labels=c(1:num_of_files), cex=0.8, col=colours)
            plot(prin_comp$x[,c(1,3)], pch='.', col=return_cutree_order[[which(n_of_groups == k_to_run)]], xlab='', ylab='')
            # text(prin_comp$x[,c(1,3)], labels=c(1:num_of_files), cex=0.8, col=return_cutree_order[[which(n_of_groups == k_to_run)]])
            text(prin_comp$x[,c(1,3)], labels=c(1:num_of_files), cex=0.8, col=colours)
            plot(prin_comp$x[,c(2,3)], pch='.', col=return_cutree_order[[which(n_of_groups == k_to_run)]], xlab='', ylab='')
            # text(prin_comp$x[,c(2,3)], labels=c(1:num_of_files), cex=0.8, col=return_cutree_order[[which(n_of_groups == k_to_run)]])
            text(prin_comp$x[,c(2,3)], labels=c(1:num_of_files), cex=0.8, col=colours)
        dev.off()

        if ( n_of_groups == 2 ){
            next;
        }
        
        if(plot_groups){
            par(mfrow = c(n_of_groups,width))
            if (width/n_of_groups > 10){
                png(paste0(experiment_name, "/", hclust_group, "/Groups_", date_time,"_", n_of_groups, "_", hclust_group, ".png"),height=width/n_of_groups*1000,width=1000,pointsize=18)
            } else {
                png(paste0(experiment_name, "/", hclust_group, "/Groups_", date_time,"_", n_of_groups, "_", hclust_group, ".png"),width=width*200,height=n_of_groups*200,pointsize=18)
            }
                
                           
                # if (width/n_of_groups > 10){
                #     mat <- t(matrix(unlist(order_of_plotting),n_of_groups,width,byrow=T))
                # } else {
                    mat <- matrix(unlist(order_of_plotting),n_of_groups,width,byrow=T)
                # }
                    
                    
                par(mar = c(rep(0.25,4)), oma = c(0,7,0,0))
                m <- layout(mat)
                layout.show(m)

                
                GroupNumberPrint <- mat[,1]
                
                for(q1 in 1:length(files)){

                    f_temp <- load_flowFrame(files=files, q1=q1, flowFrameOrPathnames = flowFrameOrPathnames)
                    plotDens(f_temp, c(xMark,yMark), cex.main=2, cex.lab=4, cex.axis=4, main="", axes=T, xlab="", ylab="", xaxt='n', yaxt='n',
                             xlim=xlim, ylim=ylim)#xlim=c(0,4), ylim=c(0,4))
                    if (length(grep(pattern = paste0("^", q1,"$"), x = GroupNumberPrint)) >= 1){
                        mtext(paste0( "Group "), x = GroupNumberPrint, side=2, line=5, cex=1,outer = F, adj=0.5)
                        mtext(paste0( grep(pattern = paste0("^", q1,"$"), x = GroupNumberPrint) ), side=2, line=1, cex=3,outer = F, adj=0.5)
                    }
                    # Add overlay 1D plot
                    if(DoOverlay){
                        overlay.dens <- density(f_temp@exprs[,xMark])
                        overlay.dens$y <- overlay.dens$y / max(overlay.dens$y)
                        if (identical(ylim, NA)) {
                            overlay.scale <- max(f_temp@exprs[,yMark])
                            overlay.translation <- min(f_temp@exprs[,yMark])
                        } else {
                            overlay.scale <- (ylim[2]-ylim[1])
                            overlay.translation <- ylim[1]
                        }
                        par(new=T)
                        lines(overlay.dens$x, overlay.dens$y*overlay.scale + overlay.translation, col="grey48")
                        
                        #horiz
                        overlay.dens <- density(f_temp@exprs[,yMark])
                        overlay.dens$y <- overlay.dens$y / max(overlay.dens$y)
                        if (identical(xlim, NA)) {
                            overlay.scale <- max(f_temp@exprs[,xMark])
                            overlay.translation <- min(f_temp@exprs[,xMark])
                        } else {
                            overlay.scale <- (xlim[2]-xlim[1])
                            overlay.translation <- xlim[1]
                        }
                        par(new=T)
                        lines(overlay.dens$y*overlay.scale + overlay.translation, overlay.dens$x, col="grey48")#)
                    }

                    # Add contour lines
                    if(DoContour){
                        z <- kde2d(f_temp@exprs[,xMark], f_temp@exprs[,yMark], n = 50)
                        # print(z)
                        z$z <- (z$z)^(1/4)
                        contour(z, drawlabels = FALSE, add = TRUE, nlevels = 10, lty = 2, col="grey49")
                    }

                    if (!is.null(vert_line)){   
                        abline(v=vert_line[q1])                    
                    }
                    if (!is.null(horiz_line)){
                        abline(h=horiz_line[q1])
                    }
                    
                   
                    # should replace with labels parameter and with upper left corner instead of relying on the xlim and ylim
                    if (!is.na(xlim) && !is.na(ylim)){
                        if (!is.null(Labels)){
                            # text(x = (xlim[1]-0.15), y = (ylim[2]-0.2), labels = Labels[q1], cex=1.75, pos = 4)
                            # text( xlim[1], ylim[2], Labels[q1], adj = c( 0.5, 1 ), col = "black") # top middle
                            text( xlim[1], ylim[2], Labels[q1] , adj = c( 0, 1 ), col = "black") # top left
                        } else{
                            text(x = (xlim[1]-0.15), y = (ylim[2]-0.2), labels = q1, cex=2, pos = 4)
                        }
                    }


                    if (!is.null(boundaries)){
                        lines(boundaries[[q1]]@boundaries, lwd=2)
                    }

                    box(lwd=2)
                }

            dev.off()
        }

        if(verbose==T){ cat("Time to plot for ", n_of_groups, " groups: ", TimeOutput(start_groups),"\n",sep="") }
    }


    names(return_cutree_order) <- names(order_of_plotting_total) <- k_to_run

    cat("Total time: ", TimeOutput(start),"\n",sep="")
    return(list(Best_k=best_k, k_scores=scores_vis, Groups=return_cutree_order, dist_all=dist_all, dist_all_norm=dist_all_norm,
                order_of_plotting=order_of_plotting_total))
}
