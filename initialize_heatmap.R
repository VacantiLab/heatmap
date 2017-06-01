
cd <- getwd()
r_scripts <- list.files(cd)
initialize_fun_index <- match('initialize_heatmap.R',r_scripts)
r_scripts <- r_scripts[-initialize_fun_index]
rank_volcano_fun_index <- match('RankVolcanoData.R',r_scripts)
r_scripts <- r_scripts[-rank_volcano_fun_index]
r_scripts_pathways <- paste(cd,r_scripts[],sep='/')
for (i in 1:length(r_scripts_pathways)) {source(r_scripts_pathways[[i]])}
