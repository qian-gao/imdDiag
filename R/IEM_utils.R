generate_diagnosis_2 <- # retired
  function( input.control = NULL, # Sample, features
            input.sample = NULL, # Sample, features
            reference = NULL, # IEM, features, no other columns
            top.nr = 3,
            calc.z = TRUE,
            top.feature.nr = 50){
    
    # Check features in new samples and reference, take all the features in reference
    f_1 <- colnames(input.sample)[!colnames(input.sample) %in% c("Sample", "IEM")]
    f_2 <- colnames(reference)[!colnames(reference) %in% c("Sample", "IEM")]
    fs <- unique(f_2)
    input.sample[, fs[!fs %in% f_1]] <- 0 
    #reference[, fs[!fs %in% f_2]] <- 0 
    
    features <- colnames(reference)[!colnames(reference) %in% c("IEM")]
    s <- input.sample[, features]
    
    # standardize
    if (calc.z) {  
      input.control[, fs[!fs %in% f_1]] <- 0
      c <- input.control[, features]
      mean.c <- apply(c, 2, mean)
      sd.c <- apply(c, 2, sd)
      
      z.s <- t(apply(s, 1, function(x){ (x - mean.c)/sd.c }))
    } else {
      z.s <- s
    }
    
    z.s[is.na(z.s)] <- 0
    reference[is.na(reference)] <- 0
    compare <-
      apply( z.s, 1,
             function(x){
               apply( reference[, features], 1,
                      function(y){
                        #cor(x,y, method = "pearson")
                        #cor(x,y, method = "spearman")
                        #(x - mean(x)) %*% (y - mean(y)) / (sqrt(sum((x - mean(x))^2))*sqrt(sum((y - mean(y))^2)))
                        index <- order(abs(t(y)), decreasing = TRUE)[1:top.feature.nr]
                        x[index] %*% y[index] / (sqrt(sum(x[index]^2))*sqrt(sum(y[index]^2)))
                        
                      })
             })
    colnames(compare) <- input.sample$Sample
    
    diagnosis <-
      data.table::rbindlist(
        lapply(1:ncol(compare), function(x){
          
          ind <- sort(compare[, x], decreasing = TRUE, index.return = TRUE)$ix[1:top.nr]
          diag <- 
            data.frame( Sample = input.sample[x, "Sample"],
                        Putative_IEM = reference$IEM[ind],
                        Score = compare[ind, x],
                        Rank = 1:top.nr)
          
          return(diag)
          
        }))
    
    diagnosis.f <-
      diagnosis %>%
      pivot_wider(names_from = "Rank", values_from = c("Putative_IEM", "Score"), names_sort = TRUE)
    
    return(diagnosis.f)
    
  }

generate_diagnosis <-
  function( input.control = NULL, # Sample, features
            input.sample = NULL, # Sample, features
            reference = NULL, # IEM, features, no other columns
            top.nr = 3,
            calc.z = TRUE){
    
    # Check features in new samples and reference, take all the features in reference
    f_1 <- colnames(input.sample)[!colnames(input.sample) %in% c("Sample", "IEM")]
    f_2 <- colnames(reference)[!colnames(reference) %in% c("Sample", "IEM")]
    # fs <- unique(c(f_1, f_2))
    # input.sample[, fs[!fs %in% f_1]] <- 0 
    # reference[, fs[!fs %in% f_2]] <- 0 
    fs <- unique(f_2)
    input.sample[, fs[!fs %in% f_1]] <- 0
    
    features <- colnames(reference)[!colnames(reference) %in% c("IEM")]
    s <- input.sample[, features]
    
    # standardize
    if (calc.z) {  
      input.control[, fs[!fs %in% f_1]] <- 0
      c <- input.control[, features]
      mean.c <- apply(c, 2, mean)
      sd.c <- apply(c, 2, sd)
      
      z.s <- t(apply(s, 1, function(x){ (x - mean.c)/sd.c }))
    } else {
      z.s <- s
    }
    
    z.s[is.na(z.s)] <- 0
    reference[is.na(reference)] <- 0
    compare <-
      apply( z.s, 1,
             function(x){
               apply( reference[, features], 1,
                      function(y){
                        #cor(x,y, method = "pearson")
                        #cor(x,y, method = "spearman")
                        #(x - mean(x)) %*% (y - mean(y)) / (sqrt(sum((x - mean(x))^2))*sqrt(sum((y - mean(y))^2)))
                        x %*% y / (sqrt(sum(x^2))*sqrt(sum(y^2)))

                      })
             })
    colnames(compare) <- input.sample$Sample
    
    diagnosis <-
      data.table::rbindlist(
        lapply(1:ncol(compare), function(x){
          
          ind <- sort(compare[, x], decreasing = TRUE, index.return = TRUE)$ix[1:top.nr]
          diag <- 
            data.frame( Sample = input.sample[x, "Sample"],
                        Putative_IEM = reference$IEM[ind],
                        Score = compare[ind, x],
                        Rank = 1:top.nr)
          
          return(diag)
          
        }))
    
    diagnosis.f <-
      diagnosis %>%
      pivot_wider(names_from = "Rank", values_from = c("Putative_IEM", "Score"), names_sort = TRUE)
    
    return(diagnosis.f)
    
  }

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

plot_iem <- function(x, iem, ref) {
  
  met <- ref[ref$IEM %in% iem, -1] 
  met <- colnames(met)[met != 0]
  
  y <- 
    x %>%
    select(Sample, IEM, all_of(met)) %>%
    filter(IEM %in% c(iem, "Control")) %>% 
    pivot_longer(-c(1:2), names_to = "Metabolite", values_to = "Intensity")
  
  
  y$IEM <- factor(y$IEM, levels = c("Control", iem))
  ggplot(y, aes(x = IEM, y = log2(Intensity))) +
    geom_boxplot(aes(color = IEM)) +
    geom_jitter(aes(color = IEM, text = Sample), alpha = 0.5) +
    facet_wrap(~ Metabolite, scales = "free_y") +
    scale_x_discrete(labels = c("Control", "Patient")) +
    scale_color_manual( values = c("#5C80BC", "#E8C547")) +
    theme_bw() +
    theme(legend.position = "top")
  
}