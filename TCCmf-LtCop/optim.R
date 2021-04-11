function (X = NULL, objfun = NULL, confun = NULL, A = NULL, B = NULL, 
          Aeq = NULL, Beq = NULL, lb = NULL, ub = NULL, tolX = 1e-05, 
          tolFun = 1e-06, tolCon = 1e-06, maxnFun = 1e+07, maxIter = 4000) 
{
  if (is.null(X)) 
    stop("please input initial value")
  if (is.null(objfun)) 
    stop("please write objective function")
  X = as.matrix(X)
  Xtarget = as.vector(X)
  dims = NULL
  dims$nxrows = nrow(X)
  dims$nxcols = ncol(X)
  dims$nvar = length(Xtarget)
  B = as.vector(B)
  Beq = as.vector(Beq)
  if (is.null(Aeq)) {
    Aeq = matrix(NA, nrow = 0, ncol = dims$nvar)
  }
  if (is.null(A)) {
    A = matrix(NA, 0, dims$nvar)
  }
  nLineareq = nrow(Aeq)
  ncolAeq = ncol(Aeq)
  nLinearIneq = nrow(A)
  Ainput = A
  lenghlb = length(lb)
  lenghub = length(ub)
  if (lenghlb > dims$nvar) {
    print("invalid bounds")
    lb = lb(1:dims$nvar)
    lenghlb = dims$nvar
  }
  else if (lenghlb < dims$nvar) {
    if (lenghlb > 0) {
      print("invalid bounds")
    }
    lb = rbind(lb, matrix(rep(-Inf, dims$nvar - lenghlb), 
                          ncol = 1))
    lenghlb = dims$nvar
  }
  if (lenghub > dims$nvar) {
    print("invalid bounds")
    ub = ub(1:dims$nvar)
    lenghub = dims$nvar
  }
  else if (lenghub < dims$nvar) {
    if (lenghub > 0) {
      print("invalid bounds")
    }
    ub = rbind(ub, matrix(rep(Inf, dims$nvar - lenghub), 
                          ncol = 1))
    lenghub = dims$nvar
  }
  len = min(lenghlb, lenghub)
  if (sum(lb > ub) > 0) {
    return("please check bounds")
  }
  start = NULL
  start$xform = X
  medx = matrix(rep(1, dims$nvar), ncol = 1)
  Xtarget[Xtarget < lb] = lb[Xtarget < lb]
  Xtarget[Xtarget > ub] = ub[Xtarget > ub]
  X = matrix(Xtarget, dim(start$xform))
  start$g = matrix(0, dims$nvar, 1)
  start$f = objfun(X)
  if (!is.null(confun)) {
    conf = confun(X)
    ctmp = conf$c
    ceqtmp = conf$ceq
    start$ncineq = as.vector(ctmp)
    start$nceq = as.vector(ceqtmp)
    start$gnc = matrix(0, dims$nvar, length(start$ncineq))
    start$gnceq = matrix(0, dims$nvar, length(start$nceq))
  }
  else {
    conf = ctmp = ceqtmp = start$ncineq = start$nceq = start$gnc = start$gnceq = NULL
  }
  if (is.null(start$ncineq)) 
    start$ncineq = matrix(0, 0, 1)
  if (is.null(start$nceq)) 
    start$nceq = matrix(0, 0, 1)
  if (is.null(start$gnc)) 
    start$gnc = matrix(0, dims$nvar, 0)
  if (is.null(start$gnceq)) 
    start.gnceq = matrix(0, dims$nvar, 0)
  fval = NULL
  lambda_out = NULL
  lambda_nc = NULL
  GRADIENT = NULL
  xform = start$xform
  iter = 0
  Xtarget = as.vector(X)
  numVar = length(Xtarget)
  DIR = matrix(1, numVar, 1)
  finalf = Inf
  steplength = 1
  HESS = diag(numVar)
  finishflag = F
  lbflag = is.finite(lb)
  ubflag = is.finite(ub)
  boundM = diag(max(lenghub, lenghlb))
  if (sum(lbflag) > 0) {
    lbM = -boundM[lbflag, 1:numVar]
    lbright = -as.matrix(lb)[lbflag, , drop = F]
  }
  else {
    lbM = NULL
    lbright = NULL
  }
  if (sum(ubflag) > 0) {
    ubM = boundM[ubflag, 1:numVar]
    ubright = as.matrix(ub)[ubflag, , drop = F]
  }
  else {
    ubM = NULL
    ubright = NULL
  }
  A = rbind(lbM, ubM, A)
  B = as.vector(c(lbright, ubright, B))
  if (length(A) == 0) {
    A = matrix(0, 0, numVar)
    B = matrix(0, 0, 1)
  }
  if (length(Aeq) == 0) {
    Aeq = matrix(0, 0, numVar)
    Beq = matrix(0, 0, 1)
  }
  LAMBDA_new = NULL
  LAMBDA = NULL
  LAMBDA_old = NULL
  X = matrix(Xtarget, dim(start$xform))
  f = start$f
  nceq = start$nceq
  ncineq = start$ncineq
  nctmp = ncineq
  nc = rbind(as.matrix(nceq), as.matrix(ncineq))
  c = rbind(rbind(rbind(Aeq %*% Xtarget - Beq, as.matrix(nceq)), 
                  as.matrix(A %*% Xtarget) - matrix(B, ncol = 1)), as.matrix(ncineq))
  nonlin_eq = length(nceq)
  nonlin_Ineq = length(ncineq)
  nLineareq = nrow(Aeq)
  nLinearIneq = nrow(A)
  eq = nonlin_eq + nLineareq
  ineq = nonlin_Ineq + nLinearIneq
  nctl = ineq + eq
  if (nonlin_eq > 0) {
    nleq_i = (1:nonlin_eq)
  }
  else {
    nleq_i = NULL
  }
  if (nonlin_Ineq > 0) {
    nlineq_i = ((nonlin_eq + 1):(nonlin_eq + nonlin_Ineq))
  }
  else {
    nlineq_i = NULL
  }
  if (eq > 0 & ineq > 0) {
    ga = rbind(abs(c[1:eq, , drop = F]), c[(eq + 1):nctl, 
                                           , drop = F])
  }
  else if (eq > 0) {
    ga = abs(c[1:eq, , drop = F])
  }
  else if (ineq > 0) {
    ga = c[(eq + 1):nctl, , drop = F]
  }
  else {
    ga = NULL
  }
  if (length(c) > 0) {
    maxga = max(ga)
  }
  else maxga = 0
  x_old = Xtarget
  c_old = c
  ggf_old = matrix(0, numVar, 1)
  ggf = start$g
  gnc = cbind(start$gnceq, start$gnc)
  tgc_old = matrix(0, nctl, numVar)
  LAMBDA = matrix(0, nctl, 1)
  lambda_nc = matrix(0, nctl, 1)
  nfval = 1
  ngval = 1
  errfloat = NULL
  while (!finishflag) {
    len_nc = length(nc)
    nctl = nLineareq + nLinearIneq + len_nc
    resfd = estgradient(Xtarget, objfun, confun, lb, ub, 
                        f, nc[nlineq_i], nc[nleq_i], 1:numVar, medx, dims, 
                        ggf, gnc[, nlineq_i, drop = F], gnc[, nleq_i, drop = F])
    ggf = resfd$gradf
    gnc[, nlineq_i] = resfd$estgIn
    gnc[, nleq_i] = resfd$estgeq
    nvals = resfd$nvals
    nfval = nfval + nvals
    if (length(gnc) > 0) {
      gc = cbind(cbind(cbind(t(Aeq), gnc[, nleq_i]), t(A)), 
                 gnc[, nlineq_i])
    }
    else if (length(Aeq) > 0 || length(A) > 0) {
      gc = cbind(t(Aeq), t(A))
    }
    else {
      gc = matrix(0, numVar, 0)
    }
    tgc = t(gc)
    if (eq > 0) {
      for (i in 1:eq) {
        iopp = tgc[i, ] %*% ggf
        if (iopp > 0) {
          tgc[i, ] = -tgc[i, ]
          c[i] = -c[i]
        }
      }
    }
    if (iter > 0) {
      maxgrad = norm(as.matrix(ggf + t(tgc) %*% lambda_nc), 
                     "I")
      if (nctl > eq) {
        maxc = norm(as.matrix(lambda_nc[(eq + 1):nctl] * 
                                c[(eq + 1):nctl]), "I")
      }
      else {
        maxc = 0
      }
      if (is.finite(maxgrad) & is.finite(maxc)) {
        errfloat = max(maxgrad, maxc)
      }
      else {
        errfloat = Inf
      }
      errga = maxga
      if (errfloat < tolFun & errga < tolCon) {
        finishflag = TRUE
      }
      else {
        if (nfval > maxnFun) {
          Xtarget = xtrial
          f = f_old
          ggf = ggf_old
          finishflag = TRUE
        }
        if (iter >= maxIter) {
          finishflag = TRUE
        }
      }
    }
    if (!finishflag) {
      iter = iter + 1
      if (ngval > 1) {
        LAMBDA_new = LAMBDA
        g_new = ggf + t(tgc) %*% LAMBDA_new
        g_old = ggf_old + t(tgc_old) %*% LAMBDA
        gdif = g_new - g_old
        xdif = Xtarget - x_old
        if (t(gdif) %*% xdif < steplength^2 * 0.001) {
          while (t(gdif) %*% xdif < -1e-05) {
            gdif[order(gdif * xdif)[1], ] = gdif[order(gdif * 
                                                         xdif)[1], ]/2
          }
          if (t(gdif) %*% xdif < (.Machine$double.eps * 
                                  norm(HESS, "F"))) {
            tgccdif = t(tgc) %*% c - t(tgc_old) %*% c_old
            tgccdif = tgccdif * (xdif * tgccdif > 0) * 
              (gdif * xdif <= .Machine$double.eps)
            weight = 0.01
            if (max(abs(tgccdif)) == 0) {
              tgccdif = 1e-05 * sign(xdif)
            }
            while (t(gdif) %*% xdif < (.Machine$double.eps * 
                                       norm(HESS, "F")) & (weight < 1/.Machine$double.eps)) {
              gdif = gdif + weight * tgccdif
              weight = weight * 2
            }
          }
        }
        if (t(gdif) %*% xdif > .Machine$double.eps) {
          HESS = HESS + (gdif %*% t(gdif))/c(t(gdif) %*% 
                                               xdif) - ((HESS %*% xdif) %*% (t(xdif) %*% 
                                                                               t(HESS)))/c(t(xdif) %*% HESS %*% xdif)
        }
      }
      else {
        LAMBDA_old = matrix(.Machine$double.eps + t(ggf) %*% 
                              ggf, nctl, 1)/(apply(t(tgc) * t(tgc), 2, sum) + 
                                               .Machine$double.eps)
        iact = 1:eq
      }
      ngval = ngval + 1
      L_old = LAMBDA
      tgc_old = tgc
      ggf_old = ggf
      c_old = c
      f_old = f
      x_old = Xtarget
      xint = matrix(0, numVar, 1)
      HESS = (HESS + t(HESS)) * 0.5
      tryCatch({
        resqp = solve.QP(HESS, -ggf, -t(tgc), c, meq = eq)
        DIR = resqp$solution
        lambda = resqp$Lagrangian
        iact = resqp$iact
      }, error = function(e) {
        resqp = solqp(HESS, ggf, tgc, -c, xint, eq, nrow(tgc), 
                      numVar)
        DIR = resqp$X
        lambda = resqp$lambda
        iact = resqp$indxact
      })
      lambda_nc[, 1] = 0
      lambda_nc[iact, ] = lambda[iact]
      lambda[(1:eq)] = abs(lambda[(1:eq)])
      if (eq > 0 & ineq > 0) {
        ga = rbind(abs(c[1:eq, , drop = F]), c[(eq + 
                                                  1):nctl, , drop = F])
      }
      else if (eq > 0) {
        ga = abs(c[1:eq, , drop = F])
      }
      else if (ineq > 0) {
        ga = c[(eq + 1):nctl, , drop = F]
      }
      if (length(c) > 0) {
        maxga = max(ga)
      }
      else maxga = 0
      LAMBDA = lambda[(1:nctl)]
      LAMBDA_old = apply(cbind(LAMBDA, 0.5 * (LAMBDA + 
                                                LAMBDA_old)), 1, max)
      ggfDIR = t(ggf) %*% DIR
      xtrial = Xtarget
      commerit = f + sum(LAMBDA_old * (ga > 0) * ga) + 
        1e-30
      if (maxga > 0) {
        commerit2 = maxga
      }
      else if (f >= 0) {
        commerit2 = -1/(f + 1)
      }
      else {
        commerit2 = 0
      }
      if (f < 0) {
        commerit2 = commerit2 + f - 1
      }
      if ((maxga < .Machine$double.eps) & (f < finalf)) {
        finalf = f
        finalx = Xtarget
        finalHess = HESS
        finalgrad = ggf
        finallambda = lambda
        finalmaxga = maxga
        finalerrfloat = errfloat
      }
      search = T
      alpha = 2
      while (search & nfval < maxnFun) {
        alpha = alpha/2
        if (alpha < 1e-04) {
          alpha = -alpha
        }
        if ((norm(as.matrix(DIR), "I") < 2 * tolX || 
             abs(alpha * ggfDIR) < tolFun) & (maxga < tolCon)) {
          finishflag = T
        }
        Xtarget = xtrial + alpha * DIR
        X = matrix(Xtarget, dim(start$xform))
        f = objfun(X)
        if (!is.null(confun)) {
          conf = confun(X)
          nctmp = conf$c
          nceqtmp = conf$ceq
        }
        else {
          conf = nctmp = nceqtmp = NULL
        }
        nctmp = as.vector(nctmp)
        nceqtmp = as.vector(nceqtmp)
        nfval = nfval + 1
        nc = c(nceqtmp, nctmp)
        c = matrix(c(Aeq %*% Xtarget - Beq, nceqtmp, 
                     A %*% Xtarget - matrix(B, ncol = 1), nctmp), 
                   ncol = 1)
        if (eq > 0 & ineq > 0) {
          ga = rbind(abs(c[1:eq, , drop = F]), c[(eq + 
                                                    1):nctl, , drop = F])
        }
        else if (eq > 0) {
          ga = abs(c[1:eq, , drop = F])
        }
        else if (ineq > 0) {
          ga = c[(eq + 1):nctl, , drop = F]
        }
        if (length(c) > 0) {
          maxga = max(ga)
        }
        else maxga = 0
        merit = f + sum(LAMBDA_old * (ga > 0) * ga)
        if (maxga > 0) {
          merit2 = maxga
        }
        else if (f >= 0) {
          merit2 = -1/(f + 1)
        }
        else {
          merit2 = 0
        }
        if (f < 0) {
          merit2 = merit2 + f - 1
        }
        search = (merit2 > commerit2) & (merit > commerit)
      }
      steplength = alpha
      if (!finishflag) {
        mf = abs(steplength)
        LAMBDA = mf * LAMBDA + (1 - mf) * L_old
        X = matrix(Xtarget, dim(start$xform))
      }
    }
  }
  GRADIENT = ggf
  if (f > finalf) {
    Xtarget = finalx
    f = finalf
    HESS = finalHess
    GRADIENT = finalgrad
    lambda = finallambda
    maxga = finalmaxga
    ggf = finalgrad
    errfloat = finalerrfloat
  }
  fval = f
  X = matrix(Xtarget, dim(start$xform))
  nLinearIneq = nrow(Ainput)
  lambda_out = NULL
  lambda_out$lower = matrix(0, lenghlb, 1)
  lambda_out$upper = matrix(0, lenghub, 1)
  if (nLineareq > 0) {
    lambda_out$eqlin = lambda_nc[1:nLineareq]
  }
  ii = nLineareq
  if (nonlin_eq > 0) {
    lambda_out$eqnonlin = lambda_nc[(ii + 1):(ii + nonlin_eq)]
  }
  ii = ii + nonlin_eq
  if (sum(lbflag != 0) > 0) {
    lambda_out$lower[lbflag] = lambda_nc[(ii + 1):(ii + sum(lbflag != 
                                                              0))]
  }
  ii = ii + sum(lbflag != 0)
  if (sum(ubflag != 0) > 0) {
    lambda_out$upper[ubflag] = lambda_nc[(ii + 1):(ii + sum(ubflag != 
                                                              0))]
  }
  ii = ii + sum(ubflag != 0)
  if (nLinearIneq > 0) {
    lambda_out$ineqlin = lambda_nc[(ii + 1):(ii + nLinearIneq)]
  }
  ii = ii + nLinearIneq
  if (nonlin_Ineq > 0) {
    lambda_out$ineqnonlin = lambda_nc[(ii + 1):length(lambda_nc)]
  }
  return(list(par = X, fn = fval, counts = cbind(nfval, ngval), 
              lambda = lambda_out, grad = GRADIENT, hessian = HESS))
}