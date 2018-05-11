######################
#
# GB 20131106
# Adapted PerformanceAnalytics code to make some changes
#
######################

require(PerformanceAnalytics)

chart.Correlation.GB = function (R, histogram = TRUE, method = c("pearson", "kendall", 
    "spearman"), ...) 
{
    x = checkData(R, method = "matrix")
    if (missing(method)) 
        method = method[1]
    panel.cor <- function(x, y, digits = 2, prefix = "", use = "pairwise.complete.obs", 
        method, cex.cor, ...) {
        usr <- par("usr")
        on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        nComplete = sum(complete.cases(x,y))
        r <- cor(x, y, use = use, method = method)
        txt <- format(c(r, 0.123456789), digits = digits)[1]
        txt <- paste(prefix, txt, sep = "")
        if (missing(cex.cor)) 
            cex.cor <- 0.5/strwidth(txt)
        test <- cor.test(x, y)
        Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
            cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", 
                "**", "*", ".", " "))
        text(0.5, 0.5, txt, cex = ifelse(cex.cor*(abs(r)) < cex.cor/10, cex.cor/10, cex.cor * abs(r)))
        text(0.8, 0.8, Signif, cex = cex.cor, col = 2)
        text(0.5, 0.2, paste("(n=", nComplete, ")", sep=""), cex = 2, col = gray(0.5))
    }
    f <- function(t) {
        dnorm(t, mean = mean(x), sd = sd.xts(x))
    }
    hist.panel = function(x, ...) {
        par(new = TRUE)
        hist(x, col = "light gray", probability = TRUE, axes = FALSE, 
            main = "", breaks = "FD")
        lines(density(x, na.rm = TRUE), col = "red", lwd = 1)
        rug(x)
    }
    if (histogram) 
        pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor, 
            diag.panel = hist.panel, method = method, ...)
    else pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor, 
        method = method, ...)
}
