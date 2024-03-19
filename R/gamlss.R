library(gamlss)

data("abdom")
abd0 <- gamlss(y~poly(x,3), data=abdom, family=NO)
summary(abd0)
predict(abd0,what = "mu",type = "response")[1]
plot(abd0)
wp(abd0)