# Load the required libraries
pacman::p_load(
  "tidyverse",
  "showtext",
  "hexSticker",
  install = FALSE
)

## Read in the dag image
img_dag <- magick::image_read(
  paste(here::here(), "inst/logo/hex-logo-dag-trimmed.png", sep = "/"),
  depth = 16
)

## Render the sticker logo
sticker(
  subplot = img_dag,
  package = "IPWBayes",
  p_size = 32,
  p_y = 1.5,
  p_family = "serif",
  p_fontface = "bold",
  s_x = 1,
  s_y = 0.93,
  s_width = 1.85,
  s_height = 1.8,
  h_fill = "#000000",
  h_color = "#8A2BE2",
  dpi = 600,
  filename="inst/logo/package-hex.png"
)
