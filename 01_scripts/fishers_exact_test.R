# Conduct Fisher's Exact Test on Chromosome Composition Data

# input: output of the plot_mod_chr_comp.R text file

# rm(list=ls())

# Set working directory
setwd("~/Documents/10_bernatchez/10_paralogs")

# Set sex of coexpression modules
#sex <- "female"
sex <- "male"

# Set filenames
in.filename <- paste(sex, c("count_of_genes_per_chr_per_module.txt"), sep = "_")

# Import
results <- read.delim2(file = in.filename, header = T, sep = "\t", row.names = 1)
#colnames(results)[1] <- "chr"

results <- results[,-c(which(colnames(results)=="baseline"))]

head(results)

#### 1. extract ####
chr.of.interest <- "NC_027304.1"
mod.of.interest <- "steelblue"

num.genes.in.chr.in.mod <- results[chr.of.interest , mod.of.interest]
num.genes.not.chr.in.mod <- sum(results[, mod.of.interest]) - num.genes.in.chr.in.mod

num.genes.in.chr.not.mod <- sum(results[chr.of.interest, ]) - num.genes.in.chr.in.mod
num.genes.not.chr.not.mod <- sum(results) - sum(results[chr.of.interest, ]) - (sum(results[, mod.of.interest]) - num.genes.in.chr.in.mod)


#### 2. stat test ####
stat.test.df <- matrix(c(num.genes.in.chr.in.mod, num.genes.not.chr.in.mod
                       , num.genes.in.chr.not.mod, num.genes.not.chr.not.mod)
                       , nrow = 2)

stat.test.df

stat.test.result <- fisher.test(stat.test.df)
str(stat.test.result)
stat.test.result$p.value


#### 3. Create loop #####
# set nulls
chr.of.interest <- NULL; mod.of.interest <- NULL; num.genes.in.chr.in.mod <- NULL; num.genes.not.chr.in.mod <- NULL
num.genes.in.chr.not.mod <- NULL; num.genes.not.chr.not.mod <- NULL; stat.test.df <- NULL
result.temp <- NULL; coll.fisher.res.temp <- NULL; coll.fisher.res.all <- NULL

for(c in rownames(results)){
  #print(c)
  
    for(m in colnames(results)){
    #print(c(c, m))
      chr.of.interest <- c
      mod.of.interest <- m
      
      #print(c(c,m))
      
      num.genes.in.chr.in.mod <- results[chr.of.interest , mod.of.interest]
      num.genes.not.chr.in.mod <- sum(results[, mod.of.interest]) - num.genes.in.chr.in.mod
      
      num.genes.in.chr.not.mod <- sum(results[chr.of.interest, ]) - num.genes.in.chr.in.mod
      num.genes.not.chr.not.mod <- sum(results) - sum(results[chr.of.interest, ]) - (sum(results[, mod.of.interest]) - num.genes.in.chr.in.mod)
      
      stat.test.df <- matrix(c(num.genes.in.chr.in.mod, num.genes.not.chr.in.mod
                               , num.genes.in.chr.not.mod, num.genes.not.chr.not.mod)
                             , nrow = 2)
      
      result.temp <- fisher.test(stat.test.df)
      coll.fisher.res.temp <- c(c, m, result.temp$p.value, num.genes.in.chr.in.mod, num.genes.not.chr.in.mod)
      
      coll.fisher.res.all <- rbind(coll.fisher.res.all, coll.fisher.res.temp)
      
      
    }
}

colnames(coll.fisher.res.all) <- c("chr", "mod", "pval", "g.in.chr.in.mod", "g.not.in.chr.in.mod")
rownames(coll.fisher.res.all) <- 1:nrow(coll.fisher.res.all)

output <- as.data.frame(coll.fisher.res.all)
output$pval <- as.numeric(as.character(output$pval))

output

# sort by ascending p-value
output <- output[with(output, order(output$pval)), ]

head(output)


#### 04. Save values data
# Write out results
filename.output <- paste(sex, c("fisher_test_chr_enrich.txt"), sep = "_")
write.table(x = output
            , file = filename.output
            , quote = F
            , sep = "\t"
            , row.names = F
            #, col.names = NA
            )


# still to do: multiple test correct

