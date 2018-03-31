## Return a design matrix at quartile (or supplied) knots, and all knots.
## The knots are taken at the quartiles, and min and max of x.
## x: all the explanatory data
## K: number of inner knots, or a vector of all knots (including extrema)
## deg: degree of the spline polynomial
TpfDesign <- function(x, K, deg) {

    res <- list()

    if (length(K) == 1 && is.numeric(K)) {
        knots <- quantile(unique(x), seq(0, 1, len = K + 2))[-c(1, K + 2)]
        names(knots) <- NULL
        res$knots <- c(min(x), knots, max(x))
    } else if (is.vector(K) && is.numeric(K)) {
        knots <- K[-c(1, length(K))]
        res$knots <- K
    } else {
        stop("Supplied knots must a numeric vector.")
    }
    names(knots) <- NULL

    splines <- outer(x, knots, `-`)
    if (deg > 0 && deg < 3) {
        splines <- splines^deg * (splines > 0)
        res$design <- cbind(1, poly(x, degree = deg, raw = T, simple = TRUE),
                            splines, deparse.level = 0)
        colnames(res$design) <- NULL
    } else {
        stop("Invalid degree.")
    }

    ## if (deg == 1) {
    ##     splines <- splines * (splines > 0)
    ##     res$design <- cbind(1, x, splines, deparse.level = 0)
    ## } else if (deg == 2) {
    ##     splines <- splines^2 * (splines > 0)
    ##     res$design <- cbind(1, x, x^2, splines, deparse.level = 0)
    ## } else {
    ##     stop("Invalid degree.")
    ## }

    return(res)
}

## Return the constraint matrix A, where Ax > 0.
## knots: all knots, i.e. extrama of x and inner knots. If deg = 1, this can
##        the number of inner knots
## shape: increasing, decreasing
## deg: degree of the polynomial splines
TpfConstMat <- function(knots, shape, deg) {
    if (deg == 1) {
        if (length(knots) == 1) {
            knots <- rep(NA, knots + 2) # a dummy knots vector
        }
        if (shape == "increasing") {
            A <- diag(length(knots) - 1)
            A[lower.tri(A)] <- 1
            A <- cbind(0, A)
            return(A)
        } else if (shape == "decreasing") {
            warning("decresing cases untested.")
            A <- diag(length(knots) - 1)
            A[lower.tri(A)] <- 1
            A <- cbind(0, A)
            return(-1 * A)
        } else {
            stop("Unrecognised shape.")
        }
    } else if (deg == 2) {
        if (!is.vector(knots, "numeric") || length(knots) == 1) {
            stop("Knots vector must be supplied for quadratic spline.")
        } else if (shape == "increasing") {
            A <- 2 * outer(knots[-(1:2)], knots[-c(1, length(knots))], `-`)
            A[upper.tri(A)] <- 0
            A <- rbind(0, 0, A, deparse.level = 0)
            A <- cbind(0, 1, 2 * knots, A, deparse.level = 0)
            return(A)
        } else if (shape == "decreasing") {
            warning("decresing cases untested.")
            A <- 2 * outer(knots[-(1:2)], knots[-c(1, length(knots))], `-`)
            A[upper.tri(A)] <- 0
            A <- rbind(0, 0, A, deparse.level = 0)
            A <- cbind(0, 1, 2 * knots, A, deparse.level = 0)
            return(-1 * A)
        } else {
            stop("Unrecognised shape.")
        }
    } else {
        stop("Unrecognised spline degree.")
    }
}

## Truncated power function, model with only a single population curve
## data : 1st col x, 2nd col y, 3rd col groups
## K : number of quantile (inner) knots, or a vector of inner knots
## deg : degree of spline polynomial
SubjectsTpf <- function(data, K, deg = 1, shape = "increasing", size = 100,
                        burn = size / 10) {

    if (deg != 1 && deg != 2) {
        stop("Invalid spline degree. Must be 1 or 2.")
    }
    if (burn < 0) {
        stop("Negative burn.")
    } else if (burn == 0) {
        warning("Not burning any samples.")
    }

    x <- data[[1]]
    y <- data[[2]]
    grp <- data[[3]]

    n.terms <- K + deg + 1
    n.samples <- NROW(data)

    ## convert the group variable into a factor
    if (is.factor(data[[3]])) {
        grp <- data[[3]]
    } else {
        grp <- factor(data[[3]], levels = unique(data[[3]]))
    }

    lvl.sub <- levels(grp)
    idx.sub <- tapply(seq_len(n.samples), grp, function(x) x)
    n.subs <- length(idx.sub)

    ## construct the design matrix and knots
    design.ls <- TpfDesign(x, K, deg)
    knots <- design.ls$knots            # all knots (with extrema)
    n.spline <- length(knots) - 2       # number of inner knots (w/o extrema)
    X.pop <- design.ls$design

    ## get rid of the design.ls to save space
    rm(design.ls)

    ## construct cross-products of the design matrix
    X.pop.sq <- crossprod(X.pop)
    X.sub.sq <- array(NA, c(n.terms, n.terms, n.subs),
                      list(NULL, NULL, lvl.sub))
    for (j in lvl.sub) {
        X.sub.sq[, , j] <- crossprod(X.pop[idx.sub[[j]], ])
    }

    ## An idempotent matrix to extract the spline terms
    Kmat <- diag(c(rep(0, deg + 1), rep(1, n.spline)))
    idx.poly <- seq_len(deg + 1)        # index of polynomial terms

    ## Construct the constraint matrix and its cross-product
    A <- TpfConstMat(knots, shape, deg)
    A.t <- t(A)

    ## Calculate an inverse of A to easily produce feasible states
    A.inv <- diag(NCOL(A))
    A.inv[row(A.inv) > diff(dim(A))] <- A
    A.inv <- solve(A.inv)

    ## initialise population coefs and subjects deviations
    kcoef.pop <- seq(0, 1, len = n.terms)
    kcoef.sub <- matrix(0, n.terms, n.subs, dimnames = list(NULL, lvl.sub))

    ## initialise prediction contribution by population coefs
    kpred.pop <- X.pop %*% kcoef.pop

    ## initialise prediction contribution by subjects deviations
    kpred.sub <- rep(NA, n.samples)
    for (j in lvl.sub) {
        idx <- idx.sub[[j]]
        kpred.sub[idx] <- X.pop[idx, ] %*% kcoef.sub[, j]
    }

    ## remove the dummy variables used in the for loop
    rm(j, idx)

    ## initialise the output list, by the order of lvl.sub
    samples <- list(population = matrix(NA, n.terms, size),
                    subjects = array(NA, c(n.terms, n.subs, size),
                                     dimnames = list(NULL, lvl.sub)),
                    precision = list())
    ## samples$precision: precisions (matrix) (num vec/mat list)
    ## pop: population coefs; poly: subject polynomial terms;
    ## sub: subject deviations; eps: error terms.
    samples$precision$pop <- rep(NA, size)
    samples$precision$poly <- array(NA, c(deg + 1, deg + 1, size))
    samples$precision$sub <- rep(NA, size)
    samples$precision$eps <- rep(NA, size)

    ## burnin followed by actual sampling
    for (k in seq.int(-burn + 1, size)) {

##### BELOW ARE THE NEW CODE THAT PRODUCES THE STRANGE BEHAVIOUR #####

        ## get the precisions (varariance-covariance matrices)
        kprecs <- TpfGetCov(list(kcoef.pop), list(kcoef.sub), list(kpred.pop),
                            list(kpred.sub), list(y), idx.poly, Kmat, 1,
                            n.subs, n.spline, n.samples)

        ## get the coefs and deviations
        kcoefs <- TpfGetCoefs(kcoef.sub, kpred.sub, X.pop, X.pop.sq, X.sub.sq,
                              lvl.sub, idx.sub, y, kprecs$eps, kprecs$pop,
                              kprecs$poly, kprecs$sub, A, A.t, A.inv, Kmat,
                              n.terms)

        ## for the ease of reading
        kcoef.pop <- kcoefs$coef.pop
        kcoef.sub <- kcoefs$coef.sub
        kpred.pop <- kcoefs$pred.pop
        kpred.sub <- kcoefs$pred.sub

        if (k > 0) {
            ## store the precisions
            samples$precision$pop[k] <- kprecs$pop
            samples$precision$poly[, , k] <- kprecs$poly
            samples$precision$sub[k] <- kprecs$sub
            samples$precision$eps[k] <- kprecs$eps

            ## store the coefs and deviations
            samples$population[, k] <- kcoef.pop
            samples$subjects[, , k] <- kcoef.sub
        }

##### ABOVE ARE THE NEW CODE THAT PRODUCES THE STRANGE BEHAVIOUR #####

    }

##### BELOW ARE THE OLD CODE THAT WORKS FINE #####
    ## ## hyperparameters of priors
    ## ig.a <- 0.001
    ## ig.b <- 0.001
    ## wi.df <- 1
    ## wi.sig <- diag(0.001, deg + 1)

    ## ## Burnin step
    ## for (k in seq_len(burn + size)) {
    ##     ## Update variances
    ##     c.vpop <- 1 / rgamma(1, shape = n.spline / 2 + ig.a,
    ##                          scale = 0.5 * crossprod(Kmat %*% c.pop) + ig.b)
    ##     c.vsub <- 1 / rgamma(1, shape = (n.subs * n.spline) / 2 + ig.a,
    ##                          scale = 0.5 * crossprod(c(Kmat %*% c.sub)) + ig.b)
    ##     c.ppoly <- rWishart(1, df = wi.df + n.subs,
    ##                         solve(wi.sig + tcrossprod(c.sub[idx.poly, ])))
    ##     resid.vec <- y - pred.pop - pred.sub
    ##     c.veps <- 1 / rgamma(1, shape = 0.5 * n.samples + ig.a,
    ##                          scale = 0.5 * crossprod(resid.vec) + ig.b)


    ##     ## Update population estimates
    ##     M.pop <- solve(X.pop.sq + (c.veps / c.vpop) * Kmat)
    ##     mu.pop <- tcrossprod(M.pop, X.pop) %*% (y - pred.sub)
    ##     sig.pop <- c.veps * M.pop

    ##     lower.pop <- A %*% c.sub
    ##     lower.pop <- -1 * apply(lower.pop, 1, min)
    ##     lower.pop[lower.pop < 0] <- 0

    ##     ## Initialise the starting values of the truncated normal sampler
    ##     start.pop <- A.inv %*% c(mu.pop[1], (lower.pop + 0.1))

    ##     c.pop <- TruncatedNormal::rmvtgauss.lin(1, mu.pop, sig.pop,
    ##                                             Amat = A.t,
    ##                                             Avec = lower.pop,
    ##                                             start = start.pop,
    ##                                             burnin = 1000)
    ##     if (k > burn) {
    ##         rec.idx <- k - burn
    ##         if (rec.idx %% 1000 == 0) {
    ##             cat(rec.idx, "samples generated. \n")
    ##         }
    ##         samples$population[, rec.idx] <- c.pop
    ##         samples$var.pop[rec.idx] <- c.vpop
    ##         samples$prec.poly[, , rec.idx] <- c.ppoly
    ##         samples$var.sub[rec.idx] <- c.vsub
    ##         samples$var.eps[rec.idx] <- c.veps
    ##     }

    ##     ## Update prediction contribution by the population curve.
    ##     pred.pop <- X.pop %*% c.pop
    ##     y.diff.pop <- y - pred.pop

    ##     ## Update subject specific estimates
    ##     lower.sub <- -1 * (A %*% c.pop)

    ##     ## Initialise the starting values of the truncated normal sampler
    ##     ## start.sub <- rep(0, n.terms)
    ##     zeros <- rep(0, n.terms)

    ##     ## Calculate the precision matrix term
    ##     G.term <- diag(1 / c.vsub, n.terms)
    ##     G.term[1:NROW(c.ppoly), 1:NCOL(c.ppoly)] <- c.ppoly
    ##     G.term <- c.veps * G.term

    ##     for (i in seq_len(n.subs)) {
    ##         idx <- idx.sub[[i]]
    ##         M.sub <- solve(X.sub.sq[[i]] + G.term)
    ##         mu.sub <- tcrossprod(M.sub, X.pop[idx, ]) %*% y.diff.pop[idx]
    ##         sig.sub <- c.veps * M.sub

    ##         c.sub[, i] <- TruncatedNormal::rmvtgauss.lin(1, mu.sub, sig.sub,
    ##                                                      Amat = A.t,
    ##                                                      Avec = lower.sub,
    ##                                                      start = zeros,
    ##                                                      burnin = 1000)

    ##         ## Update prediction contribution by subject curves
    ##         pred.sub[idx] <- X.pop[idx, ] %*% c.sub[, i]

    ##         if (k > burn) {
    ##             samples[[i]][, k - burn] <- c.sub[, i]
    ##         }
    ##     }
    ## }
##### ABOVE ARE THE OLD CODE THAT WORKS FINE #####

    ## return posterior means (coefs only), samples, and information regarding
    ## the basis functions.
    means <- list(population = rowMeans(samples$population),
                  subjects = rowMeans(samples$subjects, dims = 2))
    basis <- list(type = "tpf", knots = knots, degree = deg)
    info <- list(lvl_pop = NULL, lvl_sub = lvl.sub, n_terms = n.terms)
    data <- data.frame(x = x, y = y, grp_sub = grp, grp_pop = NA)

    list(means = means, samples = samples, basis = basis, info = info,
         data = data)
}


## Subjects model with multiple population curves
## data: 1st predictor, 2nd response, 3rd subjects, 4th populations
## K: number of quantile (inner) knots, or a vector of inner knots
## deg: degree of spline polynomial
SubjectsTpfMul <- function(data, K, deg = 1, shape = "increasing", size = 100,
                           burn = size / 10, verbose = FALSE) {
    ## if (deg != 1 && deg != 2) {
    ##     stop("Invalid spline degree. Must be 1 or 2.")
    ## }

    n.terms <- K + deg + 1              # number of coefs to describe a curve
    n.samples <- NROW(data)             # number of provided datapoints

    x <- data[[1]]
    y <- data[[2]]

    ## convert groups variables to factors
    if (is.factor(data[[3]])) {
        grp.sub <- data[[3]]
    } else {
        grp.sub <- factor(data[[3]], levels = unique(data[[3]]))
    }
    if (is.factor(data[[4]])) {
        grp.pop <- data[[4]]
    } else {
        grp.pop <- factor(data[[4]], levels = unique(data[[4]]))
    }

    ## calculate the population indices
    ## these are used to extract from the input (main) dataframe
    lvl.pop <- levels(grp.pop)
    idx.pop <- tapply(seq_len(n.samples), grp.pop, function(x) x)
    n.pops <- length(idx.pop)

    ## calculate the subject indices within each population
    ## these are used to extract from each population dataframe (P_i)
    lvl.sub <- list()
    idx.sub <- list()
    for (i in lvl.pop) {
        grp.sub.i <- grp.sub[grp.pop == i, drop = TRUE]
        lvl.sub[[i]] <- levels(grp.sub.i)
        idx.sub[[i]] <- tapply(seq_along(grp.sub.i), grp.sub.i , function(x) x)
    }
    n.subs <- vapply(idx.sub, length, 0L) # number of subjects in each pop

    ## remove the dummy variables used in the for loop
    rm(i, grp.sub.i)

    ## construct design matrix and knots
    design.ls <- TpfDesign(x, K, deg)
    knots <- design.ls$knots            # all knots (with extrema)
    n.spline <- length(knots) - 2       # number of inner knots (w/o extrema)

    ## construct cross-products of the design matrix
    ## also, seperate the design matrix for each population (P_i)
    ## X.pop.sq = crossprod(P_i), X.sub.sq = crossprod(X_ij)
    ## X.pop.sq is a 3D array, X.sub.sq is a list of 3D arrays.
    X.pop <- list()
    X.pop.sq <- list()
    X.sub.sq <- list()
    for (i in lvl.pop) {
        X.pop[[i]] <- design.ls$design[idx.pop[[i]], ]
        X.pop.sq[[i]] <- crossprod(X.pop[[i]])
        X.sub.sq[[i]] <- array(NA, c(n.terms, n.terms, n.subs[[i]]),
                               list(NULL, NULL, lvl.sub[[i]]))
        for (j in lvl.sub[[i]]) {
            X.sub.sq[[i]][, , j] <- crossprod(X.pop[[i]][idx.sub[[i]][[j]], ])
        }
    }

    ## get rid of the design.ls to save space
    rm(design.ls)

    ## get y for each population (y_i)
    y.pop <- tapply(y, grp.pop, function(x) x)

    ## an idempotent matrix to extract the spline terms
    Kmat <- diag(c(rep(0, deg + 1), rep(1, n.spline)))
    idx.poly <- seq_len(deg + 1)        # index of polynomial terms

    ## construct the constraint matrix and its cross-product
    A <- TpfConstMat(knots, shape, deg)
    A.t <- t(A)

    ## calculate an inverse of A to easily produce feasible states
    A.inv <- diag(NCOL(A))
    A.inv[row(A.inv) > diff(dim(A))] <- A
    A.inv <- solve(A.inv)

    ## initialise the current estimates
    kcoef.pop <- list()
    kcoef.sub <- list()
    kpred.pop <- list()
    kpred.sub <- list()

    for (i in lvl.pop) {
        ## initialise population coefs and subjects deviations
        kcoef.pop[[i]] <- rep(0, n.terms)
        kcoef.sub[[i]] <- matrix(0, n.terms, n.subs[[i]],
                                 dimnames = list(NULL, lvl.sub[[i]]))

        ## initialise prediction contribution by population coefs
        kpred.pop[[i]] <- X.pop[[i]] %*% kcoef.pop[[i]]

        ## initialise prediction contribution by subjects deviations
        kpred.sub[[i]] <- rep(NA, length(idx.pop[[i]]))
        for (j in lvl.sub[[i]]) {
            idx <- idx.sub[[i]][[j]]
            kpred.sub[[i]][idx] <- X.pop[[i]][idx, ] %*% kcoef.sub[[i]][, j]
        }
    }

    ## remove the dummy variables used in the for loop
    rm(i, j, idx)

    ## initialise the output list, by the order of lvl.pop and lvl.sub
    samples <- list(population = list(), subjects = list(), precision = list())
    for (i in lvl.pop) {
        ## samples$population: population coefs (num mat list)
        ## samples$subjects: subject deviations (num 3D ary list)
        samples$population[[i]] <- matrix(NA, n.terms, size)
        samples$subjects[[i]] <- array(NA, c(n.terms, n.subs[[i]], size),
                                       dimnames = list(NULL, lvl.sub[[i]]))
    }
    ## samples$precision: precisions (matrix) (num vec/mat list)
    ## pop: population coefs; poly: subject polynomial terms;
    ## sub: subject deviations; eps: error terms.
    samples$precision$pop <- rep(NA, size)
    samples$precision$poly <- array(NA, c(deg + 1, deg + 1, size))
    samples$precision$sub <- rep(NA, size)
    samples$precision$eps <- rep(NA, size)

    ## remove the dummy variables used in the for loop
    rm(i)

    ## burnin followed by actual sampling
    for (k in seq.int(-burn + 1, size)) {
        ## get the precisions (varariance-covariance matrices)
        kprecs <- TpfGetCov(kcoef.pop, kcoef.sub, kpred.pop, kpred.sub, y.pop,
                            idx.poly, Kmat, n.pops, n.subs, n.spline, n.samples)

        if (k > 0) {
            ## store the precisions
            samples$precision$pop[k] <- kprecs$pop
            samples$precision$poly[, , k] <- kprecs$poly
            samples$precision$sub[k] <- kprecs$sub
            samples$precision$eps[k] <- kprecs$eps
        }

        ## get the coefs and deviations
        kcoefs <- parallel::mcmapply(TpfGetCoefs, kcoef.sub, kpred.sub, X.pop,
                                     X.pop.sq, X.sub.sq, lvl.sub, idx.sub,
                                     y.pop,
                                     MoreArgs = list(kprecs$eps, kprecs$pop,
                                                     kprecs$poly, kprecs$sub, A,
                                                     A.t, A.inv, Kmat, n.terms),
                                     SIMPLIFY = FALSE)

        for (i in lvl.pop) {
            ## convert mapply output into a suitable format
            kcoef.pop[[i]] <- kcoefs[[i]]$coef.pop
            kcoef.sub[[i]] <- kcoefs[[i]]$coef.sub
            kpred.pop[[i]] <- kcoefs[[i]]$pred.pop
            kpred.sub[[i]] <- kcoefs[[i]]$pred.sub

            if (k > 0) {
                ## store the coefs and deviations
                samples$population[[i]][, k] <- kcoef.pop[[i]]
                samples$subjects[[i]][, , k] <- kcoef.sub[[i]]
            }
        }

        if (verbose && (k %% 1000 == 0)) {
            cat(k, " samples generated.\n")
        }
    }

    ## return posterior means (coefs only), samples, and information regarding
    ## the basis functions.
    means <- list(population = lapply(samples$population, rowMeans),
                  subjects = lapply(samples$subjects,
                                    function(x) rowMeans(x, dims = 2)))
    basis <- list(type = "tpf", knots = knots, degree = deg)
    info <- list(lvl_pop = lvl.pop, lvl_sub = lvl.sub, n_terms = n.terms)
    data <- data.frame(x = x, y = y, grp_sub = grp.sub, grp_pop = grp.pop)

    list(means = means, samples = samples, basis = basis, info = info,
         data = data)

}


## get a sample from the coefs posterior (assuming 1 population)

## different in each ITERATION and POPULATION
## coef.sub: individual deviations (num mat)
## pred.sub: prediction contribution by the sub deviations (num vec)

## different in each POPULATION
## X.pop: design matrix (num mat)
## X.pop.sq: crossprod of the whole design matrix (num mat)
## X.sub.sq: crossprod of the individual design matrix (num 3d ary)
## lvl.sub: names of each subject (str vec)
## idx.sub: indices of each subject in X.pop (num vec list)
## y.pop: response data (num vec list)

## different in each ITERATAION
## prc.eps: precision of the Gaussian noise (num)
## prc.pop: precision of the population spline terms (num)
## prc.poly: precision of the individual polynomial term (num mat)
## prc.sub: precision of the individual spline terms (num)

## constants
## A: constraint matrix (num mat)
## A.t: t(A), for the purpose of using Berwin's sampler (num mat)
## A.inv: pseudo-inverse of A for generating feasible starting values (num mat)
## Kmat: to extract spline terms (num mat)
## n.terms: number of parameters to descrive a curve (num)

## RETURN
## coef.pop: population coefs (num vec)
## coef.sub: individual deviations (num mat)
## pred.pop: prediction contribution by the pop coefs (num vec)
## pred.sub: prediction contribution by the sub deviations (num vec)
TpfGetCoefs <- function(coef.sub, pred.sub, X.pop, X.pop.sq, X.sub.sq, lvl.sub,
                        idx.sub, y.pop,
                        prc.eps, prc.pop, prc.poly, prc.sub, A,
                        A.t, A.inv, Kmat, n.terms) {

    M.pop <- solve(X.pop.sq + (prc.pop / prc.eps) * Kmat)
    mu.pop <- tcrossprod(M.pop, X.pop) %*% (y.pop - pred.sub)
    sig.pop <- M.pop / prc.eps

    lower.pop <- A %*% coef.sub
    lower.pop <- -1 * apply(lower.pop, 1, min)
    lower.pop[lower.pop < 0] <- 0

    ## Initialise the starting values of the truncated normal sampler
    start.pop <- A.inv %*% c(mu.pop[1], (lower.pop + 0.1))

    coef.pop <- TruncatedNormal::rmvtgauss.lin(1, mu.pop, sig.pop,
                                            Amat = A.t,
                                            Avec = lower.pop,
                                            start = start.pop,
                                            burnin = 1000)

    ## Update prediction contribution by the population curve.
    pred.pop <- X.pop %*% coef.pop
    y.diff.pop <- y.pop - pred.pop

    ## Update subject specific estimates
    lower.sub <- -1 * (A %*% coef.pop)

    ## initialise the starting values for the truncated normal sampler, and
    ## the coef.sub matrix
    zeros <- rep(0, n.terms)
    coef.sub <- matrix(NA, n.terms, length(lvl.sub),
                       dimnames = list(NULL, lvl.sub))


    ## Calculate the precision matrix term
    half.N <- diag(prc.sub, n.terms)
    half.N[1:NROW(prc.poly), 1:NCOL(prc.poly)] <- prc.poly
    half.N <- half.N / prc.eps

    for (j in lvl.sub) {
        idx <- idx.sub[[j]]
        M.sub <- solve(X.sub.sq[, , j] + half.N)
        mu.sub <- tcrossprod(M.sub, X.pop[idx, ]) %*% y.diff.pop[idx]
        sig.sub <- M.sub / prc.eps

        coef.sub[, j] <- TruncatedNormal::rmvtgauss.lin(1, mu.sub, sig.sub,
                                                        Amat = A.t,
                                                        Avec = lower.sub,
                                                        start = zeros,
                                                        burnin = 1000)

        ## Update prediction contribution by subject curves
        pred.sub[idx] <- X.pop[idx, ] %*% coef.sub[, j]
    }

    list(coef.pop = coef.pop, coef.sub = coef.sub, pred.pop = pred.pop,
         pred.sub = pred.sub)
}

## get a sample from the covariance/precision posterior (assuming multiple population)

## different in each ITERATION and POPULATION
## coef.pop: population coefs (num vec list)
## coef.sub: individual deviations (num mat list)
## pred.pop: prediction contribution by the pop coefs (num vec list)
## pred.sub: prediction contribution by the sub deviations (num vec list)

## different in each POPULATION
## y.pop: response data (num vec list)

## constants
## idx.poly: indices to extract polynomial terms (num vec)
## Kmat: to extract spline terms (num mat)
## n.pops: number of populations (num)
## n.subs: number of subjects (num vec)
## n.spline: number of spline terms/knots (num)
## n.samples: total sample size (num)

## RETURN
## prc.eps: precision of the Gaussian noise (num)
## prc.pop: precision of the population spline terms (num)
## prc.poly: precision of the individual polynomial term (num mat)
## prc.sub: precision of the individual spline terms (num)
TpfGetCov <- function(coef.pop, coef.sub, pred.pop, pred.sub, y.pop,
                      idx.poly, Kmat, n.pops, n.subs, n.spline, n.samples) {

    ## hyperparameters of priors
    ig.a <- 0.001
    ig.b <- 0.001
    wi.df <- 1
    wi.sig <- diag(0.001, length(idx.poly))

    if (typeof(coef.pop) == "list") coef.pop else list(coef.pop)
    if (typeof(coef.sub) == "list") coef.sub else list(coef.sub)
    if (typeof(pred.pop) == "list") pred.pop else list(pred.pop)
    if (typeof(pred.sub) == "list") pred.sub else list(pred.sub)

    ## crossprod of the spline terms
    xspl.pop <- sum(vapply(coef.pop, function(x) crossprod(Kmat %*% x), 0))
    xspl.sub <- sum(vapply(coef.sub, function(x) crossprod(c(Kmat %*% x)), 0))

    ## precision of population spline terms
    shp.pop <- (n.pops * n.spline) / 2 + ig.a
    scl.pop <- 0.5 * xspl.pop + ig.b
    prc.pop <- rgamma(1, shape = shp.pop, scale = scl.pop)

    ## precision of individual spline terms
    shp.sub <- (sum(n.subs) * n.spline) / 2 + ig.a
    scl.sub <- 0.5 * xspl.sub + ig.b
    prc.sub <- rgamma(1, shape = shp.sub, scale = scl.sub)

    ## precision of residuals
    resid.vec <- mapply(function(x, y, z) x - y - z,
                        y.pop, pred.pop, pred.sub,
                        SIMPLIFY = FALSE)
    shp.eps <- 0.5 * n.samples + ig.a
    scl.eps <- 0.5 * crossprod(unlist(resid.vec)) + ig.b
    prc.eps <- rgamma(1, shape = shp.eps, scale = scl.eps)

    ## tcrossprod of the polynomial terms, presented as a 3D array
    txpoly <- vapply(coef.sub, function(x) tcrossprod(x[idx.poly, ]),
                     diag(as.double(idx.poly)))

    ## precision of individual polynomial terms
    df.poly <- wi.df + sum(n.subs)
    Sigma.poly <- solve(wi.sig + rowSums(txpoly, dims = 2))
    prc.poly <- rWishart(1, df = df.poly, Sigma = Sigma.poly)[, , 1]

    list(pop = prc.pop, sub = prc.sub, poly = prc.poly, eps = prc.eps)
}
