
source("getMDI1.R")

dtt1<-function (phy, data, index = c("avg.sq", "avg.manhattan", "num.states"), 
    mdi.range = c(0, 1), nsim = 0, CI = 0.95, plot = TRUE, calculateMDIp = F, Ylim=c(0,2)) 
{
    disp = match.arg(index, c("avg.sq", "avg.manhattan", "num.states"))
    td <- treedata(phy, data)
    dtt.data <- geiger:::.dtt(td$phy, td$data, disp = disp)
    ltt <- sort(branching.times(td$phy), decreasing = TRUE)
    ltt <- c(0, (max(ltt) - ltt)/max(ltt))
    s <- ratematrix(td$phy, td$data)
    dtt.sims = NULL
    MDI = NULL
    ylim = c(range(pretty(dtt.data)))
    if (is.numeric(nsim)) {
        if (nsim > 0) {
            sims <- sim.char(td$phy, s, nsim)
            dtt.sims <- geiger:::.dtt(td$phy, sims)
            mean.sims <- apply(dtt.sims, 1, mean)
            median.sims <- apply(dtt.sims, 1, median)
            MDI <- unname(geiger:::.area.between.curves(ltt, apply(dtt.sims, 
                1, median), dtt.data, sort(mdi.range)))
            names(MDI) = disp
            colnames(dtt.sims) = NULL
            yy = range(dtt.sims)
            ylim = range(c(ylim, yy))
        }
    }
    if (plot) {
    	ylim=Ylim
        plot(ltt, dtt.data, xlab = "relative time", ylab = "disparity", 
            ylim = ylim, bty = "n", type = "n")
        if (!is.null(dtt.sims)) {
            poly = geiger:::.dtt.polygon(dtt.sims, ltt, alpha = 1 - CI)
            #polygon(poly[, "x"], poly[, "y"], col = geiger:::.transparency("darkgrey", 0.5), border = NA)
            polygon(poly[, "x"], poly[, "y"], col = "grey60", border = NA)

            lines(ltt, median.sims, lty = 2)
        }
        lines(ltt, dtt.data, type = "l", lwd = 1)
    }
 res = list(dtt = dtt.data, times = ltt, sim = dtt.sims, MDI = MDI)
    drp = sapply(res, function(x) is.null(x))
    if (any(drp)) 
        res = res[-which(drp)]
    if (calculateMDIp) {
        pVal <- getMDIp1(res)
        res <- c(res, MDIpVal = pVal)
    }
    return(res)
}
