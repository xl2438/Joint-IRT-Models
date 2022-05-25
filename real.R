dat_score <- read.csv("response_grade_8_final.csv")
dat_cal <- read.csv("cal_grade_8_final.csv")

head(dat_score)
head(dat_cal)

mat_score <- as.matrix(dat_score[, -1])
head(mat_score)
mat_cal <- as.matrix(dat_cal[, -1])

source("em.R")
fit <- joint_model_em(mat_score, mat_cal, 25)
