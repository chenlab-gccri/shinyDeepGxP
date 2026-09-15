# # Set the conda Python path for reticulate
# reticulate::use_python("/root/.local/share/r-miniconda/envs/r-reticulate/bin/python", required = TRUE)

#options(shiny.trace = TRUE)
#reticulate::use_condaenv("r-reticulate", required = TRUE)
reticulate::use_condaenv("r-DL_env", required = TRUE)
reticulate::py_config()