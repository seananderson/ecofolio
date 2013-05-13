### R code from vignette source 'ecofolio.Rnw'

###################################################
### code chunk number 1: ecofolio.Rnw:35-38
###################################################
library(ecofolio)
data(pinkbr)
head(pinkbr)


###################################################
### code chunk number 2: ecofolio.Rnw:44-48
###################################################
library(reshape)
library(ggplot2)
x_long <- melt(pinkbr, id.vars = "year", variable_name = "river")
ggplot(x_long, aes(year, value, colour = river)) + geom_line()


###################################################
### code chunk number 3: ecofolio.Rnw:54-55
###################################################
fit_taylor(pinkbr[,-1])


###################################################
### code chunk number 4: ecofolio.Rnw:87-88
###################################################
plot_mv(pinkbr[,-1], show = "linear", ci = TRUE)


###################################################
### code chunk number 5: ecofolio.Rnw:110-111
###################################################
pe_mv(pinkbr[,-1], ci = TRUE)


###################################################
### code chunk number 6: ecofolio.Rnw:125-126
###################################################
pe_avg_cv(pinkbr[,-1], ci = TRUE, boot_reps = 500)


###################################################
### code chunk number 7: ecofolio.Rnw:142-145
###################################################
pe_mv(pinkbr[,-1], type = "linear_robust")
pe_mv(pinkbr[,-1], type = "quadratic")
pe_mv(pinkbr[,-1], type = "linear_quad_avg")


###################################################
### code chunk number 8: ecofolio.Rnw:153-158
###################################################
par(mfrow = c(1, 2))
plot_mv(pinkbr[,-1], show = "quadratic", add_z = FALSE)
mtext("Quadratic")
plot_mv(pinkbr[,-1], show = "robust", add_z = FALSE)
mtext("Robust linear")


###################################################
### code chunk number 9: ecofolio.Rnw:167-170
###################################################
pe_mv(pinkbr[,-1], type = "linear")
pe_mv(pinkbr[,-1], type = "linear_detrended")
pe_mv(pinkbr[,-1], type = "loess_detrended")


