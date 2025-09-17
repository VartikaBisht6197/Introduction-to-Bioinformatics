# Let’s first read the raw RNAseq data in R 
RNAseq_raw = read.table("data/rawCountExpressionMatrix.txt", header = TRUE, row.names = 1,check.names=TRUE)

# Use the following command to see a small portion of the dataset
# Can you guess what 1:10 does?
# How about RNAseq_raw[1:10,1:10] ? Can you write this in “English” language?
head(RNAseq_raw[1:10,1:10])

# Open : https://github.com/VartikaBisht6197/Introduction-to-Bioinformatics/blob/main/data/Phenotypes.tsv before you procced.

# Now, let’s explore rawCountExpressionMatrix.txt or RNAseq_raw 😊
# What are the column names?
# What are the row names?
# Hm, this does not look very pretty. Why is there an “X” in front of all rownames? 
# Let’s change something in the following code to make it better. Can you guess what would it be? Use this link to figure it out 😊
# https://www.rdocumentation.org/packages/utils/versions/3.6.2/topics/read.table
RNAseq_raw = read.table("data/rawCountExpressionMatrix.txt", header = TRUE, row.names = 1,check.names=FALSE)

# Let’s see if that changed anything  
head(RNAseq_raw[1:10,1:10])
# Much better 😉

# Can you tell how many rows and columns does your dataframe RNAseq_raw  have?
dim(RNAseq_raw)
# How many genes did we measure? Ans : 61,552
# How many samples do we have? Ans : 294

# What is the entry for Gene ENSG00000225972 and sample ID 11757849? 
# Use the logic discussed previously : RNAseq_raw[1:10,1:10]
# Write the logic in “English” :
# Write the logic in R :
RNAseq_raw["",""]

# So, What is the entry for Gene ENSG00000225972 and sample ID 11757849? Ans: 9
# What does 9 mean ?

# Now that we know a little bit about how to use R and also about our dataset, lets go towards “Analysis”
# First step towards analysis is to load our DESeq2 library ( What is DESeq2 anyway ? )
library(DESeq2)

# Let’s check if we actually downloaded the DESeq2 library
# We use https://www.rdocumentation.org/packages/utils/versions/3.6.2/topics/sessionInfo to get information about all packages loaded in the session 
sessionInfo()
# Woah, woah, woah ... I do not remember downloading all these libraries ?? Huh?
# Could you guess what’s happening?
# Do you find DESeq2 ?
# What version do you have ?

# Woohoo! Now that we have answer to existential questions DESeq2 presents, let’s get it started!
# Oh wait, before just running the following commands , I want to know what kind of data does DESeq2 expect ? I mean, is it okay to just put in my raw data? Or do I have to do something ?
# Use this link to figure this out : https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# Hint : Look for information on “Input data”
...

# Let’s run the following command now :
dds = DESeqDataSetFromMatrix(countData = RNAseq_raw,
                              colData = clinical,
                              design = ~ condition)

# Welcome to the realm of errors 😊
# Let’s see, 
# Error in h(simpleError(msg, call)) : 
#  error in evaluating the argument 'x' in selecting a method for function 'nrow': object 'clinical' not found
# Hmm, makes sense, I did not define ‘clinical’ .. How would R know? 
# Let’s tell R what clinical is.. can you try to do that?
# The file is “data/clinical.txt”
# Careful with the parameters.. they might changes based on what you want to do 😊
clinical = read.table(..)

# Let’s see how it looks
head(clinical)

# Nice, now we are all set!
# Let’s run this again
dds = DESeqDataSetFromMatrix(countData = RNAseq_raw,
                              colData = clinical,
                              design = ~ condition)

# Teehee :P
# Didn’t work. I know. Let’s breath in and breath out. 
# What do you think is the problem ? Hint : https://www.biostars.org/p/9465688/ 

# So what do we see?
# Well, the rowname are the names of the sample and the column name is the “condition” we use.
# Let’s do that then.
# Let’s write it in English first :
# Then, let’s do it in R :

# What are the rownames like?
rownames(clinical)

# What do we want them to be ?
clinical$sampleID

# So, let’s assign it
rownames(clinical) = clinical$sampleID
# Can you guess why we did not use “==”

# Let’s see what it did
head(clinical)

# Okay now we only want the “Ever_smoke” column not the “sampleID” column
clinical$sampleID = NULL
# Can you try to say the logic in “English” ?
# Can you guess what happened ?

# Next, what is the last thing which we must change?
colnames(clinical) = c("condition")

# Let’s check 
head(clinical)

# So, let’s try once more
dds = DESeqDataSetFromMatrix(countData = RNAseq_raw,
                              colData = clinical,
                              design = ~ condition)
# Third time is the charm!!!
# However, there is small text there:
#  the design formula contains one or more numeric variables with integer values,
# specifying a model with increasing fold change for higher values.
# did you mean for this to be a factor? if so, first convert
# this variable to a factor using the factor() function

# Can you guess what does it mean ? hint : https://r4ds.had.co.nz/factors.html#:~:text=In%20R%2C%20factors%20are%20used,to%20work%20with%20than%20characters. 
# Let’s see what we start with
clinical$condition
# Now, lets perform the factor conversion
clinical$condition =  factor(clinical$condition)
# How is it different?
clinical$condition

# Now, let’s run the DESeq2
dds = DESeqDataSetFromMatrix(countData = RNAseq_raw,
                              colData = clinical,
                              design = ~ condition)
# No erros, No messages ... nothing .. Good job!

# Let’s check our variables
ls()
# What do you see?
# Do you recognize declaring all these variables?

# Nice, what next then? We are close to asking questions now! Can someone guess what did design = ~ condition do ?
# Let’s see how does dds look
dds

# What you are looking at is a DESeq2 object. We are yet to run DESeq2 (normalisation) function.
# But first, let’s remove genes which are lowly expressed because they would affect our normalisation.
keep <- rowMedians(counts(dds)) >= 10 

# How many are we keeping and how many are we throwing away ?
table(keep)

# Could you try to write the above command in “English” ? :
dds <- dds[keep,]

# Now that we have subsetted to only keep expressing genes let’s normalise it using DESeq
dds <- DESeq(dds)

# Ooo.. something has started..
# Can you guess what is it? ( it is not straight forward – skip it if you like 😊 )
# Hint: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts

dds 
# What are the new things which have been added to dds ?

# Now, perform the Differential gene expression test between different groups
# That is people who have smoked in their life vs people who have never smoked
# As we have only one condition to contrast with, DESeq2 assumes that is what we are asking about!
smoke_vs_nosmoke <- results(dds)

# Let's see what happened
smoke_vs_nosmoke
# What do you see?
# What do the first 2 line mean?
# What do you see in the table?
# Is the variable "smoke_vs_nosmoke" a dataframe ?


# Well, then let's change it to one 😊 
smoke_vs_nosmoke <- as.data.frame(smoke_vs_nosmoke)

# Let's check if something changed
head(smoke_vs_nosmoke)
# Where do you think the 2 line on the top go ?


# Now that we have this table, let's get plotting!!
volcano = ggplot(smoke_vs_nosmoke, aes(x = log2FoldChange, y = pvalue)) + 
            geom_point()

# Whoops!
# Error in ggplot(smoke_vs_nosmoke, aes(x = log2FoldChange, y = pvalue)) :
#  could not find function "ggplot"
# Hmm, so R does not know ggplot, let's load it then 😊 
library(ggplot2)

# Let's do this again!
volcano = ggplot(smoke_vs_nosmoke, aes(x = log2FoldChange, y = pvalue)) + 
            geom_point()

# Woohoo! Good job! But well, can't see it.
# Let's save this as a PDF then! ( do you have to do that alwaays ??)
pdf(file = "volcano.plot.pdf")
volcano
dev.off()

# Hm, does not look right. Could you guess why?


# Well, p values are generally viewed in -log10 scale 😊
volcano = ggplot(smoke_vs_nosmoke, aes(x = log2FoldChange, y = -log10(pvalue))) + 
            geom_point()
pdf(file = "volcano.plot.pdf")
volcano
dev.off()

# Nice, Now it looks like something we know !!
# Let's make a volcano plot. Why is it called a volcano ?


# Let's try to interpret the plot! Which side is smoker and which side is non smoker?
# you remember the line : 
# smoke_vs_nosmoke
# log2 fold change (MLE): condition 1 vs 0
# Wald test p-value: condition 1 vs 0
# A positive log2FC indicates that the gene or feature is expressed at a higher level in condition 1 compared to condition 0.


# Can you get the name of the first 5 gene with highest log2FoldChange ? 
smoke_vs_nosmoke <- smoke_vs_nosmoke[order(smoke_vs_nosmoke$log2FoldChange),]

# Can you get the name of the first 5 gene with lowest p vale ?
smoke_vs_nosmoke <- smoke_vs_nosmoke[order(smoke_vs_nosmoke$pvalue),]

# Is there a gene which is common in these list? ENSG00000260197 and ENSG00000229236
# Can you guess which of the points on the volcano plot are these?


# Here is a clever way to do it! Can you complete the following ?
smoke_vs_nosmoke_selected = smoke_vs_nosmoke[c("",""),]
# Let's check it
smoke_vs_nosmoke_selected

# Now , ggplot magic
volcano = ggplot(smoke_vs_nosmoke, aes(x = log2FoldChange, y = -log10(pvalue))) + 
            geom_point() + geom_point(data=smoke_vs_nosmoke_selected,aes(x = log2FoldChange, y = -log10(pvalue)), color = "red")
pdf(file = "volcano.plot.pdf")
volcano
dev.off()

# What do you see?
# Were you right?
# Can you tell the logic in "English" ?

# What are these 2 genes?






