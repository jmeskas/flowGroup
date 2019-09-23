# Written by Justin Meskas
# Last updated on September 2019
# alpha version 0.99.5

#################################################################################
# Prints out the time since start_time. Used for optimizing code.
TimeOutput <- function(start_time) {
    start_time <- as.POSIXct(start_time)
    dt <- difftime(Sys.time(), start_time, units="secs")
    # Since you only want the H:M:S, we can ignore the date...
    # but you have to be careful about time-zone issues
    format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}

load_flowFrame <- function(files=files, q1=q1, flowFrameOrPathnames = flowFrameOrPathnames){
    if (flowFrameOrPathnames == "flowFrame"){
        h1 <- files[[q1]]
    } else {
        if ( mode(files) == "list" ){
            h1 <- read.FCS(filename = paste0(files[[q1]]))
        } else {
            h1 <- read.FCS(filename = paste0(files[q1]))
        }
    }
    return(h1)
}


find_number <- function(y, order_of_plotting){
    location <- unlist(lapply(order_of_plotting, function(x) {
        temp <- which(x == y)
        if (length(temp) == 0){
            temp <- "nope"
        }
        return(temp)
    }))
    names(location) <- 1:length(location)
    location <- location[which(location != "nope")]
    return(location)
}
