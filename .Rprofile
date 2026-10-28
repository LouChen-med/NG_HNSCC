# Windows R igraph and the Conda C igraph library share a DLL basename.
# Initialize Python first so both libraries can coexist in this project.
local({
  if (.Platform$OS.type == "windows") {
    prefix <- normalizePath(file.path(R.home(), "..", ".."),
                            winslash = "/", mustWork = TRUE)
    if (tolower(basename(prefix)) == "ng_hnscc_env" &&
        file.exists(file.path(prefix, "python.exe"))) {
      .libPaths(R.home("library"))
      Sys.setenv(RETICULATE_PYTHON = file.path(prefix, "python.exe"),
                 RETICULATE_AUTOCONFIGURE = "FALSE")
      reticulate::py_run_string("import numpy, pandas, igraph, leidenalg")
    }
  }
})
