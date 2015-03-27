fdrtool <- function (x, statistic = c("normal", "correlation", "pvalue", 
    "studentt"), plot = TRUE, color.figure = TRUE, verbose = TRUE, 
    cutoff.method = c("fndr", "pct0", "locfdr"), pct0 = 0.75) 
{
    statistic = match.arg(statistic)
    cutoff.method = match.arg(cutoff.method)
    if (is.vector(x) == FALSE) 
        stop("input test statistics must be given as a vector!")
    if (length(x) < 200) 
        warning("There may be too few input test statistics for reliable FDR calculations!")
    if (statistic == "pvalue") {
        if (max(x) > 1 & min(x) < 0) 
            stop("input p-values must all be in the range 0 to 1!")
    }
    if (verbose) 
        cat("Step 1... determine cutoff point\n")
    if (cutoff.method == "pct0") {
        if (statistic == "pvalue") 
            x0 = quantile(x, probs = 1 - pct0)
        else x0 = quantile(abs(x), probs = pct0)
    }
    else if (cutoff.method == "locfdr" & (statistic == "normal" | 
        statistic == "correlation")) {
        if (statistic == "normal") 
            z = x
        if (statistic == "correlation") 
            z = atanh(x)
        iqr = as.double(diff(quantile(z, probs = c(0.25, 0.75))))
        sdhat = iqr/(2 * qnorm(0.75))
        N = length(z)
        b = ifelse(N > 5e+05, 1, 4.3 * exp(-0.26 * log(N, 10)))
        z0 = b * sdhat
        if (statistic == "normal") 
            x0 = z0
        if (statistic == "correlation") 
            x0 = tanh(z0)
    }
    else {
        if (cutoff.method == "locfdr") 
            warning("cutoff.method=\"locfdr\" only available for normal and correlation statistic.")
        x0 = fndr.cutoff(x, statistic)
    }
    if (verbose) 
        cat("Step 2... estimate parameters of null distribution and eta0\n")
    cf.out <- censored.fit(x = x, cutoff = x0, statistic = statistic)
    if (statistic == "pvalue") 
        scale.param = NULL
    else scale.param <- cf.out[1, 5]
    eta0 = cf.out[1, 3]
    if (verbose) 
        cat("Step 3... compute p-values and estimate empirical PDF/CDF\n")
    nm = get.nullmodel(statistic)
    pval = nm$get.pval(x, scale.param)
    ee <- ecdf.pval(pval, eta0 = eta0)
    g.pval <- grenander(ee)
    f.pval = approxfun(g.pval$x.knots, g.pval$f.knots, method = "constant", 
        rule = 2)
    f0.pval = function(x) return(ifelse(x > 1 | x < 0, 0, rep(1, 
        length(x))))
    F.pval = approxfun(g.pval$x.knots, g.pval$F.knots, method = "linear", 
        yleft = 0, yright = g.pval$F.knots[length(g.pval$F.knots)])
    F0.pval = function(x) return(ifelse(x > 1, 1, ifelse(x < 
        0, 0, x)))
    fdr.pval = function(p) pmin(eta0/f.pval(p), 1)
    Fdr.pval = function(p) pmin(eta0 * p/F.pval(p), 1)
    if (verbose) 
        cat("Step 4... compute q-values and local fdr\n")
    qval <- Fdr.pval(pval)
    lfdr <- fdr.pval(pval)
    result = list(pval = pval, qval = qval, lfdr = lfdr, statistic = statistic, 
        param = cf.out)
    if (plot) {
        if (verbose) 
            cat("Step 5... prepare for plotting\n")
        if (statistic == "pvalue") {
            f0 <- function(zeta) return(nm$f0(zeta, scale.param))
            F0 <- function(zeta) return(nm$F0(zeta, scale.param))
            get.pval <- function(zeta) return(nm$get.pval(1 - 
                zeta, scale.param))
            x0 = 1 - x0
        }
        else {
            f0 <- function(zeta) return(2 * nm$f0(zeta, scale.param))
            F0 <- function(zeta) return(2 * nm$F0(zeta, scale.param) - 
                1)
            get.pval <- function(zeta) return(nm$get.pval(zeta, 
                scale.param))
        }
        fdr = function(zeta) fdr.pval(get.pval(zeta))
        Fdr = function(zeta) Fdr.pval(get.pval(zeta))
        F = function(zeta) 1 - eta0 * get.pval(zeta)/Fdr(zeta)
        FA = function(zeta) (F(zeta) - eta0 * F0(zeta))/(1 - 
            eta0)
        f = function(zeta) eta0 * (f0(zeta))/fdr(zeta)
        fA = function(zeta) (f(zeta) - eta0 * f0(zeta))/(1 - 
            eta0)
        ax = abs(x)
        if (statistic == "pvalue") 
            ax = 1 - ax
        xxx = seq(0, max(ax), length.out = 500)
        ll = pvt.plotlabels(statistic, scale.param, eta0)
        par(mfrow = c(3, 1))
        if (color.figure) 
            cols = c(2, 4)
        else cols = c(1, 1)
        hist(ax, freq = FALSE, bre = 50, main = ll$main, xlab = ll$xlab, 
            cex.main = 1.8)
        lines(xxx, eta0 * f0(xxx), col = cols[1], lwd = 2, lty = 3)
        lines(xxx, (1 - eta0) * fA(xxx), col = cols[2], lwd = 2)
        if (statistic == "pvalue") 
            pos1 = "topleft"
        else pos1 = "topright"
        legend(pos1, c("Mixture", "Null Component", "Alternative Component"), 
            lwd = c(1, 2, 2), col = c(1, cols), lty = c(1, 3, 
                1), bty = "n", cex = 1.5)
        plot(xxx, F(xxx), lwd = 1, type = "l", ylim = c(0, 1), 
            main = "Density (first row) and Distribution Function (second row)", 
            xlab = ll$xlab, ylab = "CDF", cex.main = 1.5)
        lines(xxx, eta0 * F0(xxx), col = cols[1], lwd = 2, lty = 3)
        lines(xxx, (1 - eta0) * FA(xxx), col = cols[2], lwd = 2)
        plot(xxx, Fdr(xxx), type = "l", lwd = 2, ylim = c(0, 
            1), main = "(Local) False Discovery Rate", ylab = "Fdr and fdr", 
            xlab = ll$xlab, lty = 3, cex.main = 1.5)
        lines(xxx, fdr(xxx), lwd = 2)
        if (eta0 > 0.98) 
            pos2 = "bottomleft"
        else pos2 = "topright"
        legend(pos2, c("fdr (density-based)", "Fdr (tail area-based)"), 
            lwd = c(2, 2), lty = c(1, 3), bty = "n", cex = 1.5)
        par(mfrow = c(1, 1))
        rm(ax)
    }
    if (verbose) 
        cat("\n")
    return(result)
}
