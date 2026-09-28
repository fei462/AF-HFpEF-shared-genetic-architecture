# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP10A_FINAL_multi_layer_integration.R
# See repository README.md for execution order and external dependencies.


# ============================================================
# STEP10A_FINAL
# Multi-layer genetic evidence integration
#
# Frozen FINAL INPUT ONLY:
# D:/A/data/STEP10A/STEP10A_INPUT_FINAL
#
# Required subfolders:
#
# STEP10A_INPUT_FINAL/
# ├──01_conjFDR/
# ├──02_FUMA/
# ├──03_LAVA/
# ├──04_coloc/
# ├──05_SuSiE/
# ├──06_eQTL_SMR/
# └──07_sQTL_SMR/
#
# Output:
# D:/A/data/STEP10A/STEP10A_FINAL_RESULTS/
#
# ============================================================

rm(list=ls())
options(stringsAsFactors=FALSE, scipen=999)

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

INPUT_ROOT <- "D:/A/data/STEP10A/STEP10A_INPUT_FINAL"

OUTPUT_ROOT <- "D:/A/data/STEP10A/STEP10A_FINAL_RESULTS"

DIRS <- c(
"00_QC",
"01_source_data",
"02_main_tables",
"03_supplementary_tables",
"04_main_figures",
"05_supplementary_figures",
"06_logs"
)

for(x in DIRS){
dir.create(
file.path(OUTPUT_ROOT,x),
recursive=TRUE,
showWarnings=FALSE
)
}


# ------------------------------------------------------------
# Helper
# ------------------------------------------------------------

find_file <- function(folder, pattern){

files <- list.files(
folder,
pattern=pattern,
recursive=TRUE,
full.names=TRUE
)

if(length(files)==0){
return(NA)
}

files[1]
}


safe_read <- function(x){

if(is.na(x) || !file.exists(x)){
return(data.table())
}

fread(x)
}


# ------------------------------------------------------------
# 1. Import final datasets
# ------------------------------------------------------------


# conjFDR
conj_files <- list.files(
file.path(INPUT_ROOT,"01_conjFDR"),
pattern="conjfdr_0.05_loci.csv",
recursive=TRUE,
full.names=TRUE
)

conj <- rbindlist(
lapply(conj_files,fread),
fill=TRUE
)

# FUMA
fuma_files <- list.files(
file.path(INPUT_ROOT,"02_FUMA"),
pattern="\\.txt$|\\.csv$",
recursive=TRUE,
full.names=TRUE
)

fuma_list <- lapply(
fuma_files,
function(x){
tryCatch(
fread(x),
error=function(e)data.table()
)
}
)

fuma <- rbindlist(
fuma_list,
fill=TRUE
)


# LAVA

lava_file <- find_file(
file.path(INPUT_ROOT,"03_LAVA"),
"local_rg_all"
)

lava <- safe_read(lava_file)


# coloc

coloc_file <- find_file(
file.path(INPUT_ROOT,"04_coloc"),
"key_results"
)

coloc <- safe_read(coloc_file)


# SuSiE

susie_file <- find_file(
file.path(INPUT_ROOT,"05_SuSiE"),
"SER_.*csv|credible"
)

susie <- safe_read(susie_file)


# eQTL SMR

eqtl_file <- find_file(
file.path(INPUT_ROOT,"06_eQTL_SMR"),
"SMR_HEIDI_results"
)

eqtl <- safe_read(eqtl_file)


# sQTL SMR

sqtl_file <- find_file(
file.path(INPUT_ROOT,"07_sQTL_SMR"),
"FULLGENOME_sQTL_SMR_HEIDI_results"
)

sqtl <- safe_read(sqtl_file)



# ------------------------------------------------------------
# 2. QC summary
# ------------------------------------------------------------

qc <- data.table(
module=c(
"conjFDR",
"FUMA",
"LAVA",
"coloc",
"SuSiE",
"eQTL_SMR",
"sQTL_SMR"
),

rows=c(
nrow(conj),
nrow(fuma),
nrow(lava),
nrow(coloc),
nrow(susie),
nrow(eqtl),
nrow(sqtl)
)
)

fwrite(
qc,
file.path(
OUTPUT_ROOT,
"00_QC",
"STEP10A_FINAL_input_QC.csv"
)
)



# ------------------------------------------------------------
# 3. Regulatory evidence layer
# ------------------------------------------------------------

eqtl_primary <- eqtl[
Primary_pass==TRUE
]

sqtl_primary <- sqtl[
Primary_pass==TRUE
]


regulatory <- rbindlist(
list(

eqtl_primary[
,.(layer="eQTL_SMR",
trait,
tissue,
gene=Gene,
locus=paste0(
Probe_Chr,":",Probe_bp
),
p_value=p_SMR,
HEIDI=p_HEIDI)
],

sqtl_primary[
,.(layer="sQTL_SMR",
trait,
tissue,
gene=Gene,
locus=paste0(
Probe_Chr,":",Probe_bp
),
p_value=p_SMR,
HEIDI=p_HEIDI)
]

),
fill=TRUE
)


fwrite(
regulatory,
file.path(
OUTPUT_ROOT,
"02_main_tables",
"STEP10A_regulatory_master_evidence.csv"
)
)



# ------------------------------------------------------------
# 4. Integrated locus evidence table
# ------------------------------------------------------------

# create unified evidence matrix

gene_summary <- regulatory[
!is.na(gene),
.(

n_regulatory_layers=
uniqueN(layer),

regulatory_layers=
paste(
unique(layer),
collapse=";"

),

traits=
paste(
unique(trait),
collapse=";"

),

tissues=
paste(
unique(tissue),
collapse=";"

),

minimum_P=
min(
p_value,
na.rm=TRUE
)

),
by=gene
]


fwrite(
gene_summary,
file.path(
OUTPUT_ROOT,
"02_main_tables",
"STEP10A_final_gene_evidence_summary.csv"
)
)


# ------------------------------------------------------------
# 5. Source tables backup
# ------------------------------------------------------------

save_list <- list(
conjFDR=conj,
FUMA=fuma,
LAVA=lava,
coloc=coloc,
SuSiE=susie,
eQTL_SMR=eqtl,
sQTL_SMR=sqtl
)


for(n in names(save_list)){

fwrite(
save_list[[n]],
file.path(
OUTPUT_ROOT,
"01_source_data",
paste0(
"FINAL_",
n,
"_source.csv"
)
)
)

}



# ------------------------------------------------------------
# 6. Main Figure
# Evidence landscape
# ------------------------------------------------------------

plot_data <- gene_summary[
order(
minimum_P
)
][1:min(30,.N)]


p <- ggplot(
plot_data,
aes(
x=reorder(
gene,
-log10(minimum_P)
),
y=-log10(minimum_P),
size=n_regulatory_layers
)
)+

geom_point()+

coord_flip()+

theme_classic(
base_size=14
)+

labs(
x=NULL,
y="-log10(P)",
size="Evidence layers"
)


ggsave(
file.path(
OUTPUT_ROOT,
"04_main_figures",
"Figure_STEP10A_evidence_landscape.tiff"
),
p,
width=7,
height=6,
dpi=600,
compression="lzw"
)



# ------------------------------------------------------------
# 7. completion
# ------------------------------------------------------------

capture.output(
sessionInfo(),
file=file.path(
OUTPUT_ROOT,
"06_logs",
"sessionInfo.txt"
)
)


cat(
"
================================================
STEP10A FINAL COMPLETED

Check:
00_QC
02_main_tables
03_supplementary_tables
04_main_figures

================================================
"
)
