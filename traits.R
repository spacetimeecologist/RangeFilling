bien.height <- function(Sciname){
  bientrait <- 'Height'
  bienname <- c('TraitName', 'TraitValue', 'Unit')
  biendf <- bientraits[which(bientraits$taxonAuthorVerbatim==Sciname), match(bienname, colnames(bientraits))]
  biendf <- biendf[is.na(match(biendf$TraitName, bientrait))!=TRUE,]
  if(nrow(biendf)>1){
    biendf.mean <- biendf[1,]
    biendf.mean$TraitValue <- mean(as.numeric(biendf$TraitValue))
    biendf <- biendf.mean
  }
  height <- as.numeric(biendf$TraitValue)
  return(height)
}

bien.form <- function(Sciname){
  bientrait <- 'Growth form'
  bienname <- c('TraitName', 'TraitValue', 'Unit')
  biendf <- bientraits[which(bientraits$taxonAuthorVerbatim==Sciname), match(bienname, colnames(bientraits))]
  biendf <- biendf[is.na(match(biendf$TraitName, bientrait))!=TRUE,]
  form <- biendf$TraitValue
  return(form)
}

traits.bien <- function(Sciname, bientraits) {
  # this block of code finds Growth Form and Height from the bientraits table
  bienlist <- list('Growth form', 'Height')
  bienname <- c('TraitName', 'TraitValue', 'Unit')
  biendf <- bientraits[which(bientraits$taxonAuthorVerbatim==Sciname), match(bienname, colnames(bientraits))]
  biendf <- biendf[is.na(match(biendf$TraitName, bienlist))!=TRUE,]
  biendf.clip <- data.frame(stringsAsFactors=FALSE)
  
  # this loop calculates means of dupilicates
  if(anyDuplicated(biendf$TraitName)>0){
    dup.index <- duplicated(biendf$TraitName)
    dup.name <- as.character(biendf$TraitName[dup.index])[1]
    dup.values <- as.character(biendf$TraitValue[biendf$TraitName==dup.name])
    dup.values <- as.numeric(dup.values)
    dup.avg <- mean(dup.values)
    biendf.clip <- biendf[-c(which(dup.index)), match(bienname, colnames(biendf))]
    biendf.clip$TraitValue[biendf.clip$TraitName==dup.name] <- dup.avg
    biendf <- biendf.clip
  }
  
  biendf <- data.frame(0,2, stringsAsFactors=FALSE)
  colnames(biendf) <- list('GrowthForm', 'Height')
  biendf$Height <- biendf.clip$TraitValue[which(biendf.clip$TraitName=='Height')]
  biendf$GrowthForm <- biendf.clip$TraitValue[which(biendf.clip$TraitName=='Growth form')]
  return(biendf) 
}



traits.try <- function(Sciname, trait, trytraits) {
  # this block of code finds species traits from the trytraits table
  trydf <- trytraits[which(gsub(' ', '', x=trytraits$SpeciesName)==Sciname),]
  trydf <- trydf[which(trydf$TraitName==trait),]
  
  # checks to make sure species has traits in try table, if not sets entry blank
  if(nrow(trydf)!=0){
    
    trydf.clip <- data.frame('TraitName'=as.character(trait), 
                             'Value'='retrieving',
                             'Unit'=toString(trydf$UnitName[which(trydf$TraitName==trait)][1]),
                             stringsAsFactors=F)
    
    if(nrow(trydf)>1){
      # checks if data type is numerical or text. If text, combines entries into one string
      if(all(is.na(as.integer(trydf$OrigValueStr[which(trydf$TraitName==trait)])))==TRUE){
        trydf.clip$Value <- toString(trydf$OrigValueStr[which(trydf$TraitName==trait)])
        trydf.clip$Unit <- toString(trydf$OrigUnitStr[which(trydf$TraitName==trait)])
      }
      # else, create a mean of the numerical values
      else{
        # checks if one of the entries is a 'mean' value, if so take mean of means
        if(is.na(match(x=TRUE, table=(grepl(pattern='mean', x=
                                              trydf$OriglName[which(trydf$TraitName==trait)],
                                            ignore.case=TRUE)
        )))==FALSE){
          
          trydf.clip$Value <- mean(trydf$StdValue[grep(pattern='mean', 
                                                       x=trydf$OriglName[], 
                                                       ignore.case=TRUE)])
          trydf.clip$Unit <- trydf$UnitName[grep(pattern='mean', 
                                                 x=trydf$OriglName[], 
                                                 ignore.case=TRUE)][1]
        }
        # else, create mean
        else{
          # checks to see if standardized values are available for mean
          if(any(is.na(trydf$StdValue[which(trydf$TraitName==trait)]))==FALSE){
            trydf.clip$Value <- mean(trydf$StdValue[which(trydf$TraitName==trait)])
            trydf.clip$Unit <- trydf$UnitName[which(trydf$Trait==trait)][1]
          }
          # else, use orignal values for mean
          else{
            trydf.clip$Value <- mean(trydf$OrigValueStr[which(trydf$TraitName==trait)])
            trydf.clip$Unit <- trydf$OrigUnitStrUnitName[which(trydf$TraitName==trait)]
          }
        }
      }
    }
    
    # else, just copy single entry
    if(nrow(trydf)==1){
      # checks if standardized values are available
      if(any(is.na(trydf$StdValue[which(trydf$TraitName==trait)]))==FALSE){
        trydf.clip$Value <- toString(trydf$StdValue[which(trydf$TraitName==trait)])
        trydf.clip$Unit <- toString(trydf$UnitName[which(trydf$TraitName==trait)])
      }
      # else, use original values
      else{
        trydf.clip$Value <- toString(trydf$OrigValueStr[which(trydf$TraitName==trait)])
        trydf.clip$Unit <- trydf$OrigUnitStr[which(trydf$TraitName==trait)]
      }
    }
    
    trydf.mod <- data.frame(paste(trydf.clip$Value, trydf.clip$Unit, sep=' '),
                            stringsAsFactors=F)
    colnames(trydf.mod) <- trydf.clip$TraitName
    trydf <- trydf.mod
  }
  
  if(nrow(trydf)==0){
    trydf <- data.frame(NA)
    colnames(trydf) <- as.character(trait) 
  }
  
  return(trydf)
}
