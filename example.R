source("periodic.R")
sitka10 <- read.table("sitka10.txt", header = T)
sitka10 <- with(sitka10, data.frame(x = days / 674,
                                    y = log.size,
                                    grp.sub = id.num,
                                    grp.pop = ozone))

## Reorder the group names just for fun
sitka10$grp.sub <- factor(sitka10$grp.sub,
                          levels = c("1", "2", "3", "4", "5",
                                     "60", "59", "56", "57", "58"))

## Model with a single population curve (about 200sec to run)
system.time(fm1 <- SubjectsTpf(sitka10, 5, deg = 2, shape = "increasing", size = 100))

## Model with multiple population curves (about 120sec to run if using linux/mac)
system.time(fm2 <- SubjectsTpfMul(sitka10, 5, deg = 2, shape = "increasing", size = 100))

## Use fm1$samples and fm2$samples to examine the samples



