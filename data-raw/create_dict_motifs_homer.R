######################################
## CREATION OF UNIFIED MOTIF DATASET
######################################
##---------------------------------------------------------------
## Input: output of HOMER findMotifsGenome.pl, knownResults.txt. 
## Includes a list of 323 (RAR:RXR 2 times) motifs and its the 
## same of all the program runs.
##---------------------------------------------------------------
## Output: data.frame containing the Motif.Name altogether with
## the ensembl_id, entrez_id and gene_symbol of the TF.
##---------------------------------------------------------------

###############################################################################
### Source: knownMotifs.txt
### Found in: After running findMotifsGenome.pl; yourfolder/knownResults.txt
###----------------------------------------------------------------------------
## Loading complete list of motifs and extract gene symbols
tf <- read.delim("~/tools/maRge/data-raw/knownResults.txt", stringsAsFactors=F)
tf <- unique(tf$Motif.Name)
tf_head <- unlist(sapply(strsplit(tf, "/"), function(x) x[[1]]))
tf_gene <- unlist(sapply(strsplit(tf_head, "(", fixed=TRUE), function(x) x[[1]]))

## Obtain possible aliases
queryGeneNames <- tf_gene
dbCon <- org.Hs.eg.db::org.Hs.eg_dbconn()
# write your SQL query
sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
# execute the query on the database
aliasSymbol <- DBI::dbGetQuery(dbCon, sqlQuery)
result <- aliasSymbol[toupper(aliasSymbol$alias_symbol) %in% 
                        toupper(as.character(queryGeneNames)),c(2,5)]

## Merge into single dataframe
tf.df <- data.frame("Motif.Name"=tf,
                    "Motif.Head"=tf_head,
                    "Motif.Symbol"=tf_gene,
                    "alias_symbol"=toupper(tf_gene),
                    stringsAsFactors=FALSE)
result$alias_symbol <- toupper(result$alias_symbol)
tf.df <- merge(tf.df, result, all.x=TRUE)
tf.df <- tf.df[,-1]

###############################################################################
### Source: table.txt
### Found in: homer/motifs/extras/table.txt
###----------------------------------------------------------------------------
## Get table from homer/motifs/extras/table.txt
table <- read.delim("~/tools/maRge/data-raw/table.txt", stringsAsFactors=F)
table$Motif.Head <- paste0(table$Factor.Name, "(", table$DBD, ")")
table <- table[,c(2,13,4,12)]
colnames(table) <- c("Motif.Name", "Motif.Head", "Motif.Symbol", "symbol")



###############################################################################
### Source: all.motifs
### Found in: homer/data/knownTFs/vertebrates/all.motifs
###----------------------------------------------------------------------------
## Get JASPAR motifs
jasp <- "grep '^>.*MA.*Jaspar*' ~/tools/maRge/data-raw/all.motifs" # Obtain HS Jaspar motifs
jasp.mot <- system(jasp, intern=TRUE)
jasp.mot <- gsub(">", "", gsub("\t0", "", jasp.mot))
jasp.mot <- unlist(sapply(strsplit(jasp.mot, "\t", fixed=TRUE), function(x) x[[2]]))

jasp.id <- unlist(sapply(strsplit(jasp.mot, "/", fixed=TRUE), function(x) x[[2]]))[-1]

## Obtain symbols from JASPAR
opts <- list()
opts[["matrixtype"]] <- "PWM"
opts[["all_versions"]] <- TRUE
PFMatrixList <- TFBSTools::getMatrixSet(JASPAR2016::JASPAR2016, opts)

jasp.df <- data.frame("Motif.Name"=jasp.mot[-1],
                      "Motif.Head"=jasp.id,
                      "Motif.Symbol"=NA)

for (i in 1:nrow(jasp.df)) {
  jasp.df$Motif.Symbol[i] <- name(PFMatrixList[[grep(jasp.df$Motif.Head[i], names(PFMatrixList))]])
}

jasp.df$alias_symbol <- toupper(jasp.df$Motif.Symbol)

## Search for alias
result.jasp <- aliasSymbol[toupper(aliasSymbol$alias_symbol) %in% 
                        toupper(as.character(jasp.df$Motif.Symbol)),c(2,5)]



result.jasp$alias_symbol <- toupper(result.jasp$alias_symbol)
jasp.df <- merge(jasp.df, result.jasp, all.x=TRUE)
jasp.df <- jasp.df[,-1]

###############################################################################
### Merge all three datasets into 1
###----------------------------------------------------------------------------
tf.all <- unique(rbind(tf.df, jasp.df, table))
tf.all <- unique(tf.all[,-2])

###############################################################################
### Convert from symbol to ID
###----------------------------------------------------------------------------
ensembl = biomaRt::useEnsembl(biomart = "ensembl", GRCh = 37, 
                              dataset = "hsapiens_gene_ensembl")
dat.hs = biomaRt::getBM(attributes = c("external_gene_name",
                                       "ensembl_gene_id",
                                       "chromosome_name"), 
                        filters = "external_gene_name", 
                        values = tf.all$symbol,
                        mart = ensembl)
colnames(tf.all)[3] <- "external_gene_name" 
tf.all <- merge(tf.all, dat.hs) ## Removing all those rows that do not have ensembl_id

tf.all <- tf.all[tf.all$chromosome_name %in% c(1:22, "X", "Y"),] # Remove non-cannonical chr
tf.all <- tf.all[,c(2,3,1,4)]
dictHomerMotifs <- tf.all

devtools::use_data(dictHomerMotifs, overwrite=TRUE)