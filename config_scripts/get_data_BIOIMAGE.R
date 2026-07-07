if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)

options(timeout = 1200)

# Base URL for the BioImage Archive accession hosting the CODEX TMA dataset
bioimage_base_url <- "https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BIAD/612/S-BIAD3612/Files"

# Directories to store the downloaded CODEX data
codex_dir <- here("data", "CODEX")
mask_dir <- file.path(codex_dir, "segmentation_masks")
boundary_dir <- file.path(codex_dir, "tumor_boundary_masks")
for (d in c(codex_dir, mask_dir, boundary_dir)) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
}
if (!dir.exists(here("aux_data"))) {
  dir.create(here("aux_data"))
}


# Download Per-cell Marker Tables -------------------------------------------------------------

download.file(file.path(bioimage_base_url, "analysis_tables/Instanseg_final_high_quality_markers.csv"),
              destfile = file.path(codex_dir, "Instanseg_final.csv"), mode = "wb")
download.file(file.path(bioimage_base_url, "analysis_tables/nimbus_cell_table_high_quality_markers.csv"),
              destfile = file.path(codex_dir, "nimbus_cell_table_all.csv"), mode = "wb")


# Derive the Core List and Marker Metadata from the Archive's Metadata Tables ------------------

# The archive doesn't include standalone cores.txt/marker_info.csv files -- derive them instead
# from the deposited core_metadata.tsv and marker_panel.tsv metadata tables
core_metadata <- read.delim(file.path(bioimage_base_url, "metadata/core_metadata.tsv"))
core_ids <- sub("^Core_", "", core_metadata$core_id)
writeLines(core_ids, con = file.path(codex_dir, "cores.txt"))

# Rename the marker panel's columns to the marker_info.csv schema expected by CODEX_TMA_Final.R
marker_panel <- read.delim(file.path(bioimage_base_url, "metadata/marker_panel.tsv"))
marker_info <- data.frame(Marker = marker_panel$analysis_column_name,
                          Prim_Location = marker_panel$primary_location,
                          Second_Location = marker_panel$secondary_location,
                          Cell_Type = marker_panel$cell_type,
                          Cell_State = marker_panel$cell_state)
write.csv(marker_info, file = here("aux_data", "marker_info.csv"), row.names = FALSE, quote = FALSE)


# Download Per-core Segmentation and Tumor-boundary Annotations -------------------------------

for (core_id in core_ids) {
  message(paste0("Downloading annotations for core: ", core_id))
  download.file(file.path(bioimage_base_url, paste0("annotations/cell_nuclei_segmentation_geojson/Core_", core_id, "_cell_nuclei.geojson")),
                destfile = file.path(mask_dir, paste0("Core_", core_id, "_mask.geojson")), mode = "wb")
  download.file(file.path(bioimage_base_url, paste0("annotations/tumor_boundary_geojson/Core_", core_id, "_tumor_boundary.geojson")),
                destfile = file.path(boundary_dir, paste0("Core_", core_id, "_tumor.geojson")), mode = "wb")
}

message("CODEX data successfully downloaded and organized in the project data directory: ", codex_dir,
       "\nMarker metadata written to: ", here("aux_data/marker_info.csv"))
