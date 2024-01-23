
# In R:
      # 2- Functional Divergence mean distance from centroid
      if (!is.na(match('FDiv', FDmetric))){
        vert1 <- scan("vert.txt", quiet = T)
        vert2 <- vert1 + 1
        vertices <- vert2[-1]
        trvertices <- spectraits[vertices, ]
        # coordinates of the center of gravity of the vertices (Gv)
        baryv <- apply(trvertices, 2, mean)
        # euclidian dstances to Gv (dB)
        distbaryv <- rep(0, nbSpecies)
        for (j in 1:nbSpecies) distbaryv[j] <- (sum((spectraits[j, ] - baryv)^2) ) ^0.5
        # mean of dB values
        meandB <- mean(distbaryv)
        # deviations to mean of db
        devdB <- distbaryv - meandB
        # computation of FDiv
        FDiv <- (sum(devdB/nbSpecies) + meandB) / (sum(abs(devdB/nbSpecies)) + meandB)
      }
    }
    if (!is.na(match('FEve', FDmetric))){
      # computation of minimum spanning tree and conversion of the 'mst' matrix into 'dist' class
      tr.dist <- stats::dist(spectraits)
      linkmst <- ape::mst(tr.dist)
      mstvect <- as.dist(linkmst)
      # computation of EW for the (nbSpecies - 1) segments to link the nbSpecies points
      EW <- rep(0, nbSpecies - 1)
      flag <- 1
      for (m in 1 : ((nbSpecies - 1) * nbSpecies / 2)) {
        if (mstvect[m] != 0) {
          EW[flag] <- tr.dist[m]
          flag <- flag + 1
        }
      }
