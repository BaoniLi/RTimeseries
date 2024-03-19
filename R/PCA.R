data <- read.csv('E:/PCA/PCA_input/18_CQ1.csv', header = T)
rh.pr <- princomp(data, cor = T)
summary(rh.pr, loadings = T)
pca_data <- predict(rh.pr)
