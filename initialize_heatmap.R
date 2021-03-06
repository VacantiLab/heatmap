
cd <- getwd()
r_scripts <- list.files(cd)
initialize_fun_index <- match('initialize_heatmap.R',r_scripts)
r_scripts <- r_scripts[-initialize_fun_index]
notes_index <- match('Notes.txt',r_scripts)
r_scripts <- r_scripts[-notes_index]
python_fun_index <- match('make_gdf.py',r_scripts)
r_scripts <- r_scripts[-python_fun_index]
r_scripts_pathways <- paste(cd,r_scripts[],sep='/')
for (i in 1:length(r_scripts_pathways))
{
    #print(i)
    source(r_scripts_pathways[[i]])
}
