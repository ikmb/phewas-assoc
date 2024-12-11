library(data.table)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
regenie_sum <- fread(args[1]) %>% mutate(p.value=10**-LOG10P)
output_name <- args[2]

"estlambda" <- function(data, plot=FALSE, proportion=1.0,
                        method="regression", filter=TRUE, df=1,... ) {
        data <- data[which(!is.na(data))]
        if (proportion>1.0 || proportion<=0)
                stop("proportion argument should be greater then zero and less than or equal to one")

        ntp <- round( proportion * length(data) )
        if ( ntp<1 ) stop("no valid measurements")
        if ( ntp==1 ) {
                warning(paste("One measurement, lambda = 1 returned"))
                return(list(estimate=1.0, se=999.99))
        }
        if ( ntp<10 ) warning(paste("number of points is too small:", ntp))
        if ( min(data)<0 ) stop("data argument has values <0")
        if ( max(data)<=1 ) {
#		lt16 <- (data < 1.e-16)
#		if (any(lt16)) {
#			warning(paste("Some probabilities < 1.e-16; set to 1.e-16"))
#			data[lt16] <- 1.e-16
#		}
                data <- qchisq(data, 1, lower.tail=FALSE)
        }
        if (filter)
        {
                data[which(abs(data)<1e-8)] <- NA
        }
        data <- sort(data)
        ppoi <- ppoints(data)
        ppoi <- sort(qchisq(ppoi, df=df, lower.tail=FALSE))
        data <- data[1:ntp]
        ppoi <- ppoi[1:ntp]
#	s <- summary(lm(data~offset(ppoi)))$coeff
#       bug fix thanks to Franz Quehenberger

        out <- list()
        if (method=="regression") {
                s <- summary( lm(data~0+ppoi) )$coeff
                out$estimate <- s[1,1]
                out$se <- s[1,2]
        } else if (method=="median") {
                out$estimate <- median(data, na.rm=TRUE)/qchisq(0.5, df)
                out$se <- NA
        } else if (method=="KS") {
                limits <- c(0.5, 100)
                out$estimate <- estLambdaKS(data, limits=limits, df=df)
                if ( abs(out$estimate-limits[1])<1e-4 || abs(out$estimate-limits[2])<1e-4 )
                        warning("using method='KS' lambda too close to limits, use other method")
                out$se <- NA
        } else {
                stop("'method' should be either 'regression' or 'median'!")
        }

        if (plot) {
                lim <- c(0, max(data, ppoi,na.rm=TRUE))
#		plot(ppoi,data,xlim=lim,ylim=lim,xlab="Expected",ylab="Observed", ...)
                oldmargins <- par()$mar
                par(mar=oldmargins + 0.2)
                plot(ppoi, data,
                     xlab=expression("Expected " ~ chi^2),
                     ylab=expression("Observed " ~ chi^2),
                     ...)
                abline(a=0, b=1)
                abline(a=0, b=out$estimate, col="red")
                par(mar=oldmargins)
        }

        out
}

gg_qqplot <- function(ps, ci = 0.95, thin=T, thin.obs.places=2) {

  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )

  if(thin){
    df <- df %>% mutate(observed=observed %>% round(thin.obs.places),
                  expected=expected %>% round(thin.obs.places)) %>% unique()
  }
  
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}

inflation <- function(chisq) {
  #chisq <- qchisq(1 - ps, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}

ps <- regenie_sum$p.value[!is.na(regenie_sum$p.value)]
qq_plot <- gg_qqplot(ps, thin = T, thin.obs.places=2) +
  theme_bw(base_size = 24) +
  annotate(
    geom = "text",
    x = -Inf,
    y = Inf,
    hjust = -0.15,
    vjust = 1 + 0.15 * 3,
    label = sprintf("λ = %.2f", estlambda(regenie_sum$p.value, method = "median")$estimate),
    size = 8
  ) +
  theme(
    axis.ticks = element_line(linewidth = 0.5),
    panel.grid = element_blank()
    # panel.grid = element_line(size = 0.5, color = "grey80")
  )


ggsave(filename=paste0(output_name,"_qq.jpg"),plot=qq_plot, height=7, width = 7)