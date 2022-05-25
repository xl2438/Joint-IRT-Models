dat <- read.csv("reponse_grade_8.csv")
dat_cal <- read.csv("cal_grade_8.csv")

id_nr <- NULL
for (i in 1:dim(dat)[1]) {
  flag <- 9 %in% dat[i, ]
  if (flag) id_nr <- c(id_nr, i)
}
dat <- dat[-id_nr, ]
id_to_remove <- NULL
for (i in 1:dim(dat)[1]) {
  flag <- 7 %in% dat[i, ]
  if (flag) id_to_remove <- c(id_to_remove, i)
}
dat <- dat[-id_to_remove, ]
id_to_remove <- NULL
for (i in 1:dim(dat)[1]) {
  flag <- 8 %in% dat[i, ]
  if (flag) id_to_remove <- c(id_to_remove, i)
}
dat <- dat[-id_to_remove, ]
dim(dat)
item_poly <- c(16, 19, 20)
dat[, item_poly] <- 1 * (dat[, item_poly] != 1)
head(dat)
ind_na <- NULL
for (i in 1:dim(dat_cal)[1]) {
  flag <- sum(is.na(dat_cal[i, ]))
  if (flag != 0) ind_na <- c(ind_na, i)
}
ind_na
dat <- dat[-ind_na, ]

dat_cal <- read.csv("cal_grade_8.csv")
head(dat_cal)
dim(dat_cal)
dim(dat)
dat_cal <- dat_cal[dat_cal$BookletNumber %in% dat$BookletNumber, ]
dim(dat_cal)

ind_na <- NULL
for (i in 1:dim(dat_cal)[1]) {
  flag <- sum(is.na(dat_cal[i, ]))
  if (flag != 0) ind_na <- c(ind_na, i)
}
ind_na

mat_cal <- as.matrix(dat_cal[, -1])
apply(mat_cal, 2, sum)

write.csv(dat, "response_grade_8_final.csv", row.names = FALSE)
write.csv(dat_cal, "cal_grade_8_final.csv", row.names = FALSE)

dat_score <- read.csv("response_grade_8_final.csv")
dat_cal <- read.csv("cal_grade_8_final.csv")
dim(dat_cal)
dim(dat_score)

ind_remove <- which(apply(is.na(dat_cal), 1, sum) != 0)
ind_remove2 <- which(apply(is.na(dat_score), 1, sum) != 0)

to_remove <- c(ind_remove, ind_remove2)
length(to_remove)

dat_cal <- dat_cal[-to_remove, ]
dat_score <- dat_score[-to_remove, ]

sample_ind <- sample(c(1:dim(dat_cal)[1]), 3000)

sample_cal <- dat_cal[sample_ind, ]
sample_score <- dat_score[sample_ind, ]
write.csv(sample_cal, "sample_cal.csv", row.names = FALSE)
write.csv(sample_score, "sample_score.csv", row.names = FALSE)

sample_cal <- read.csv("sample_cal.csv")
sample_score <- read.csv("sample_score.csv")

dim(sample_cal)
dim(sample_score)

mat_cal <- as.matrix(sample_cal[, -1])
mat_score <- as.matrix(sample_score[, -1])

apply(mat_cal, 2, sum)
apply(mat_score, 2, sum)

ind_remove <- which(apply(mat_cal, 2, sum) < 50)
mat_cal <- mat_cal[, -ind_remove]
mat_score <- mat_score[, -ind_remove]

write(colnames(mat_cal), "accNum_list.csv")

w <- t(mat_cal)
x <- t(mat_score)

library(Rcpp)
sourceCpp("test.cpp")
fit <- estimate(x, w)
fit$rho
fit$beta
fit$b
fit$ga
#############################################################

################## plots ######################
v_cal <- apply(mat_cal, 2, sum)
v_score <- apply(mat_score, 2, sum)
ind <- 1:length(v_cal)
png("cal_use_by_item.png", width = 5, height = 5, units = "in", res = 500)
plot(v_cal ~ ind, xlab = "item", ylab = "number of uses")
dev.off()
###############################
png("number_correct_by_item.png", width = 5, height = 5, units = 'in', res =
 500)
plot(v_score ~ ind, xlab = "item", ylab = "number of correct answers")
dev.off()
##################################
i_cal <- apply(mat_cal, 1, sum)
i_score <- apply(mat_score, 1, sum)
i_ind <- c(1:3000)
png("cal_use_by_person.png", width = 5, height = 5, units = "in", res = 500)
boxplot(i_cal, ylab = "number of uses")
dev.off()
png("number_correct_by_person.png", width = 5, height = 5, units = "in", res =
 500)
boxplot(i_score, ylab = "number of correct answers")
dev.off()
##################################
library(MASS)
k <- kde2d(i_score, i_cal)
png("cor_density.png", width = 5, height = 5, units = "in", res = 500)
image(k, xlab = "raw sum score", ylab = "calculator use")
dev.off()
##################################
ind <- which(i_score == 10 | i_score == 9)
j <- 4
m_cal <- mat_cal[ind, j]
m_score <- mat_score[ind, j]
table(m_cal)
table(m_cal, m_score)
###############################

alpha <- fit$alpha
beta <- fit$beta
a <- fit$a
b <- fit$b
ga <- fit$ga
rho <- fit$rho

j <- 15
alpha_se <- fit$se[1:j]
beta_se <- fit$se[(j+1):(2*j)]
a_se <- fit$se[(2*j+1):(3*j)]
b_se <- fit$se[(3*j+1):(4*j)]
ga_se <- fit$se[(4*j+1):(5*j)]
rho_se <- fit$se[5*j+1]

mat_ga <- cbind(c(1:15), ga, ga_se)
colnames(mat_ga) <- c("item", "estimate", "s.e.")
library(stargazer)
stargazer(mat_ga, out = "temp.tex")

library(ltm)
fit_x <- ltm(mat_score ~ z1)
fit_w <- ltm(mat_cal ~ z1)
fit_x$log.Lik
fit_w$log.Lik
logLik_ind <- fit_x$log.Lik + fit_w$log.Lik
sta <- -2 * (logLik_ind + 44635.4)
1 - pchisq(sta, df = 15 + 1)