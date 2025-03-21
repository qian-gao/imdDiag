generate_diagnosis <-
  function( input.control = NULL, # Sample, features
            input.sample = NULL, # Sample, features
            reference = NULL, # IMD, features, no other columns
            top.nr = 3,
            calc.z = TRUE){
    
    # Check features in new samples and reference, take all the features in reference
    f_1 <- colnames(input.sample)[!colnames(input.sample) %in% c("Sample", "IMD")]
    f_2 <- colnames(reference)[!colnames(reference) %in% c("Sample", "IMD")]
    fs <- unique(f_2)
    input.sample[, fs[!fs %in% f_1]] <- 0
    
    features <- colnames(reference)[!colnames(reference) %in% c("IMD")]
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
                        Putative_IMD = reference$IMD[ind],
                        Score = compare[ind, x],
                        Rank = 1:top.nr)
          
          return(diag)
          
        }))
    
    diagnosis.f <-
      diagnosis %>%
      pivot_wider(names_from = "Rank", values_from = c("Putative_IMD", "Score"), names_sort = TRUE)
    
    return(diagnosis.f)
    
  }