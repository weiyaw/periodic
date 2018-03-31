fm <- readRDS("periodic.rds")
str(fm)

## two populations, five subjects in each population
## quadratic spline model with truncated power series bases (five knots)
## one curve for each population (modelled as fixed)
## one curve for each subject (modelled as random)

## the population curve of Population `1`
## the intecerpt
plot(fm$population$`1`[1, 1:100])

## the quadratic term
plot(fm$population$`1`[3, 1:100])

## the first spline term
plot(fm$population$`1`[4, 1:100])

## the last (fifth) spline term
plot(fm$population$`1`[8, 1:100])


## the deviation of Subject "2" in Population `1`
## the intecerpt
plot(fm$subjects$`1`[1, "2", 1:100])

## the quadratic term
plot(fm$subjects$`1`[3, "2", 1:100])

## the first spline term
plot(fm$subjects$`1`[4, "2", 1:100])

## the last (fifth) spline term
plot(fm$subjects$`1`[8, "2", 1:100])


## variance of the population spline mixed effect (\sigma_u)
plot(1 / fm$precision$pop[1:100])

## variance of the subject intercept mixed effect (\Sigma_b[1, 1])
plot(1 / fm$precision$poly[1, 1, 1:100])

## variance of the subject spline mixed effect (\sigma_v)
plot(1 / fm$precision$sub[1:100])

## variance of the noise (\sigma_\eps)
plot(1 / fm$precision$eps[1:100])
