library(cropcircles)
library(ggpath)
library(showtext)
library(ggplot2)

# --- 1. Crop Square Plot into Hexagon ---
# Takes the 4x4 square PNG and crops it into a hex frame
hex_crop(
  images = "voronoi.png",
  to = "spStackCOS_hex.png",
  border_colour = "black",
  border_size = 15,
  just = "center"
)

# --- 2. Setup Google Fonts ---
font_add_google("Ubuntu", "ubuntu")
font_add_google("Roboto", "roboto")
showtext_auto()

ft <- "ubuntu"
pkg_name <- "spStackCOS"
txt_color <- "black"  # Use "black" or "white" depending on background contrast

# --- 3. Assemble Hex Sticker with Text Overlay ---
ggplot() +
  geom_from_path(aes(0.5, 0.5, path = "spStackCOS_hex.png")) +
  annotate("text", x = 0.5, y = 0.62, label = pkg_name, family = ft, size = 58,
           fontface = "bold", colour = txt_color, hjust = 0.5, lineheight = 0.25) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void() +
  coord_fixed() +
  theme(plot.background = element_blank(), plot.margin = margin(0, 0, 0, 0))

# --- 4. Save Final Hex Logo (SVG / PNG) ---
ggsave("spStackCOS_hex.svg", bg = NULL,
       height = 2, width = 1.73, units = "in", dpi = 900)

ggsave("spStackCOS_hex.png", bg = "transparent",
       height = 2, width = 1.73, units = "in", dpi = 900)
