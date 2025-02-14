variable_selection <-
  function( input.control = NULL, # Sample, IEM, features
            input.sample = NULL, # Sample, IEM, features
            reference = NULL, # IEM, features, no other columns, standardized
            top.nr = 20,
            start.run = FALSE,
            calc.z = TRUE
  ){
    # Check features in new samples and reference, take union
    if (!is.null(reference)) {
      f_1 <- colnames(input.control)[!colnames(input.control) %in% c("Sample", "IEM")]
      f_2 <- colnames(reference)[!colnames(reference) %in% c("Sample", "IEM")]
      fs <- unique(c(f_1, f_2))
      input.control[, fs[!fs %in% f_1]] <- 0 
      input.sample[, fs[!fs %in% f_1]] <- 0 
      reference[, fs[!fs %in% f_2]] <- NA 
      
    }
    
    if (start.run) {
      features <- colnames(input.control)[!colnames(input.control) %in% c("Sample", "IEM")]
      
      if (calc.z) {
        mean.c <- apply(input.control[, features], 2, mean)
        sd.c <- apply(input.control[, features], 2, sd)
        
        z.ref <- 
          cbind(input.sample[, c("Sample", "IEM")], 
                t(apply(input.sample[, features], 1, function(x){ (x - mean.c)/sd.c })))
        
      } else {
        z.ref <- 
          cbind(input.sample[, c("Sample", "IEM")], 
                input.sample[, features])
        
      }
      
      reference <- 
        z.ref %>%
        select(-Sample)
    }
    
    features <- colnames(input.control)[!colnames(input.control) %in% c("Sample", "IEM")]
    iem <- unique(input.sample$IEM)
    
    # standardize
    c <- as.matrix(input.control[, features])
    
    if (calc.z) {
      mean.c <- apply(c, 2, mean)
      sd.c <- apply(c, 2, sd)
      
      z.c <- t(apply(c, 1, function(x){ (x - mean.c)/sd.c }))
      z.s <- t(apply(as.matrix(input.sample[, features]), 1, 
                     function(x){ (x - mean.c)/sd.c }))
    } else {
      z.c <- c
      z.s <- as.matrix(input.sample[, features])
    }
    
    z.out <- data.frame( matrix( nrow = length(iem), ncol = length(features)))
    for (i in 1:length(iem)){
      
      iem.i <- iem[i]
      subset <- z.s[ input.sample$IEM %in% iem.i, , drop = FALSE]
      ref <- reference[ !reference$IEM %in% iem.i, features]
      
      # Intersection of selected features within same IEM
      if (nrow(subset) > 1){
        
        ws <- 
          apply(subset, 1, 
                function(x){
                  perm.out <- sparcl::HierarchicalSparseCluster.permute( rbind(x, z.c), 
                                                                         wbounds=c(1.5,2:10),
                                                                         nperms=5 )
                  
                  shc <- sparcl::HierarchicalSparseCluster( dists=perm.out$dists,
                                                            wbound=perm.out$bestw, 
                                                            method="centroid" )
                  
                  cut.off <- sort(shc$ws, decreasing = TRUE)[top.nr]
                  ws <- ifelse(shc$ws > cut.off, shc$ws, 0)
                  
                  return(ws)
                })
        
      } else {
        
        perm.out <- sparcl::HierarchicalSparseCluster.permute( rbind(subset, z.c), 
                                                               wbounds=c(1.5,2:10),
                                                               nperms=5 )
        
        shc <- sparcl::HierarchicalSparseCluster( dists=perm.out$dists,
                                                  wbound=perm.out$bestw, 
                                                  method="centroid" )
        
        cut.off <- sort(shc$ws, decreasing = TRUE)[top.nr]
        ws <- ifelse(shc$ws >= cut.off, shc$ws, 0)
        
      }
      
      
      subset.m <- apply(subset, 2, mean)
      z.iem.i <- sapply(1:nrow(ws), 
                        function(x){ 
                          ifelse(sum(ws[x,] == 0) == 0, subset.m[x], 0)
                        })
      
      # Unite of selected features with other IEM
      ws.ref <- 
        apply(ref, 1, 
              function(x){
                perm.out <- sparcl::HierarchicalSparseCluster.permute( rbind(x, 
                                                                             x, 
                                                                             subset.m
                                                                             #z.c[sample(1:5, 1), ]
                ), 
                wbounds=c(1.5,2:10),
                nperms=5 )
                
                shc <- sparcl::HierarchicalSparseCluster( dists=perm.out$dists,
                                                          wbound=perm.out$bestw, 
                                                          method="centroid" )
                
                cut.off <- sort(shc$ws, decreasing = TRUE)[1]
                ws <- ifelse(shc$ws > cut.off, shc$ws, 0)
                
                return(ws)
              })   
      
      if (ncol(ws.ref) > 1){
        
        z.iem.i.ref <- sapply(1:nrow(ws.ref), 
                              function(x){ 
                                ifelse(z.iem.i[x] != 0 | max(ws.ref[x,]) > 0, subset.m[x], 0)
                              })
        
      } else {
        
        z.iem.i.ref <- z.iem.i
        
      }
      
      # Take the average between new samples and the reference
      ref.orig <- reference[ reference$IEM %in% iem.i, features]
      if (nrow(ref.orig) > 0 && !start.run){
        
        newset <- rbind(ref.orig, z.iem.i.ref)
        z.iem.i.ref.merge <-
          apply(newset, 2, function(x){ #mean(x)
            ifelse(sum(x == 0) == 0 | sum(is.na(x)) == 1,
                   mean(x, na.rm = TRUE),
                   # {x[1] <- x[1]*0.8
                   # x[-1] <- x[-1]*0.2
                   # mean(x, na.rm = TRUE)}, # 20240924 options to do when continuous new samples coming
                   0)
          })
        
        if (sum(z.iem.i.ref.merge != 0) < 20){ ## 20240924: changed 5 to 20, maybe delete the whole section
          z.iem.i.ref.merge <-
            apply(newset, 2, function(x){ #mean(x)
              ifelse(!is.na(x[1]) && x[1] != 0, 
                     mean(x, na.rm = TRUE),
                     # {x[1] <- x[1]*0.8
                     # x[-1] <- x[-1]*0.2
                     # mean(x, na.rm = TRUE)}, # 20240924 options to do when continuous new samples coming
                     0)
            })## 20240924: changed x[1] > 0 to x[1] != 0 
          
        }
        
      } else {
        
        z.iem.i.ref.merge <- z.iem.i.ref
        
      }
      
      z.out[i, ] <- z.iem.i.ref.merge
      
    }
    
    r <- cbind( IEM = iem, z.out)
    
    colnames(r) <- c("IEM", features)
    
    rest <- reference[ !reference$IEM %in% iem, ]
    
    ref.new <- rbind(rest, r)
    
    return(ref.new)
    
  }
