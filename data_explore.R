dat <- read.csv("response_grade_8.csv")

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

