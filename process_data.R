dat <- read.csv("2017_math_statetime_score_calculator.csv")
dat <- dat[, c("AccNum", "BlockCode", "BookletNumber", "CalculatorStatus",
 "Grade", "Score")]
dat <- unique(dat)
dim(dat)
dat_cal_on <- dat[dat$CalculatorStatus == "On", ]
dim(dat_cal_on)
dat <- dat[, -4]
dat <- unique(dat)
dim(dat)
dat_p <- merge(dat, dat_cal_on, by = c("AccNum", "BookletNumber", "BlockCode",
"Grade", "Score"), all.x = TRUE)
head(dat_p)
dat_p$CalculatorStatus[is.na(dat_p$CalculatorStatus)] <- "Off"
head(dat_p)
table(dat_p$CalculatorStatus)

dat_1 <- dat_p[dat_p$Grade == 4, ]
dat_2 <- dat_p[dat_p$Grade == 8, ]

head(dat_2)
dat_score <- dat_2[, c("BookletNumber", "AccNum", "Score")]
dat_cal <- dat_2[, c("BookletNumber", "AccNum", "CalculatorStatus")]
dat_score <- tidyr::spread(dat_score, AccNum, Score)
head(dat_score)
dim(dat_score)
dat_cal[, "CalculatorStatus"] <- 1 * (dat_cal[, "CalculatorStatus"] == "On")
head(dat_cal)
dat_cal <- tidyr::spread(dat_cal, AccNum, CalculatorStatus)
head(dat_cal)
dim(dat_cal)



write.csv(dat_score, "reponse_grade_8.csv", row.names=FALSE)
write.csv(dat_cal, "cal_grade_8.csv", row.names=FALSE)

apply(dat_cal[, -1], 2, sum)