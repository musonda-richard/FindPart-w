library(tidyverse)
library(igraph)
library(ggraph)
library(ggrepel)
library(grid)


bestsha <- as.character(read.csv("../01-USA-data-analysis/rgs-AIC.csv")[1, 2])

df_est <- read.csv(paste0("../01-USA-data-analysis/", bestsha, "_estimates.csv"))
df_blocks <- read.csv(paste0("../01-USA-data-analysis/", bestsha, "_blocks.csv"))

rho_map <- df_est$k
names(rho_map) <- df_est$variant

block_map <- df_blocks$block
names(block_map) <- df_blocks$variant

# UNALIAS MAP
unalias <- c(
  "JN.1"      = "B.1.1.529.2.86.1.1",
  "JN.1.11"   = "B.1.1.529.2.86.1.1.11",
  "JN.1.16"   = "B.1.1.529.2.86.1.1.16",
  "JN.1.16.1" = "B.1.1.529.2.86.1.1.16.1",
  "JN.1.18.6" = "B.1.1.529.2.86.1.1.18.6",
  "JN.1.40"   = "B.1.1.529.2.86.1.1.40",
  "JN.1.8"    = "B.1.1.529.2.86.1.1.8",
  "KP.1.1.3"  = "B.1.1.529.2.86.1.1.11.1.1.1.3",
  "KP.2"      = "B.1.1.529.2.86.1.1.11.1.2",
  "KP.2.3"    = "B.1.1.529.2.86.1.1.11.1.2.3",
  "KP.3"      = "B.1.1.529.2.86.1.1.11.1.3",
  "KP.3.1"    = "B.1.1.529.2.86.1.1.11.1.3.1",
  "KP.3.1.1"  = "B.1.1.529.2.86.1.1.11.1.3.1.1",
  "KP.3.3"    = "B.1.1.529.2.86.1.1.11.1.3.3",
  "LP.8.1"    = "B.1.1.529.2.86.1.1.11.1.1.1.3.8.1",
  "MB.1.1"    = "B.1.1.529.2.86.1.1.49.1.1.1",
  "MC.1"      = "B.1.1.529.2.86.1.1.11.1.3.1.1.1"
)

# SAVE PANGO DISTANCE MATRIX FOR REFERENCE
pango_distance <- function(a, b) {
  
  a_parts <- strsplit(unalias[a], "\\.")[[1]]
  b_parts <- strsplit(unalias[b], "\\.")[[1]]
  
  common <- 0
  
  for(i in seq_len(min(length(a_parts), length(b_parts)))) {
    if(a_parts[i] == b_parts[i]) {
      common <- common + 1
    } else {
      break
    }
  }
  
  (length(a_parts) - common) + (length(b_parts) - common)
}

lineages <- names(unalias)

dist_mat <- matrix(0, nrow = length(lineages), ncol = length(lineages), dimnames = list(lineages, lineages))

for(i in seq_along(lineages)) {
  for(j in seq_along(lineages)) {
    dist_mat[i, j] <- pango_distance(lineages[i], lineages[j])
  }
}

write.csv(as.data.frame(dist_mat), "pango_distance_matrix.csv", row.names = TRUE)

# RECOMBINANT LINKS
recomb_edges <- tribble(
  ~from, ~to,
  "KP.3.3",   "XEC",
  "KS.1.1",   "XEC",
  "KP.3.2",   "XDY",
  "LB.1.2.1", "XDY"
)

dummy_edges <- tribble(
  ~from, ~to,
  "Dummy2_XEC", "Dummy1_XEC",
  "Dummy1_XEC", "KS.1.1",
  "Dummy2_XDY", "Dummy1_XDY",
  "Dummy1_XDY", "KP.3.2",
  "Dummy1_XDY", "LB.1.2.1"
)

# FUNCTION TO GENERATE ANCESTORS FROM JN.1
get_ancestors <- function(x) {
  
  parts <- strsplit(x, "\\.")[[1]]
  out <- vector()
  
  for(i in seq_along(parts)) {
    out[i] <- paste(parts[1:i], collapse = ".")
  }
  
  root <- "B.1.1.529.2.86.1.1"
  root_index <- which(out == root)
  
  out[root_index:length(out)]
}


# BUILD STRICT PANGO TREE
edges_full <- list()

for(v in names(unalias)) {
  
  anc <- get_ancestors(unalias[v])
  
  if(length(anc) >= 2) {
    for(i in 2:length(anc)) {
      edges_full[[length(edges_full) + 1]] <- tibble(from = anc[i - 1], to = anc[i])
    }
  }
}

edges_full <- bind_rows(edges_full) %>% distinct()


# NODE TABLE

all_nodes <- unique(c(edges_full$from, edges_full$to))

nodes <- tibble(name = all_nodes)

nodes$observed <- nodes$name %in% unalias

reverse_map <- names(unalias)
names(reverse_map) <- unalias

nodes$display_name <- reverse_map[nodes$name]
nodes$rho <- rho_map[nodes$display_name]
nodes$block <- block_map[nodes$display_name]

nodes$dummy <- !nodes$observed
nodes$recombinant <- FALSE
nodes$ghost <- FALSE

nodes$label <- ifelse(
  nodes$dummy,
  "",
  ifelse(
    is.na(nodes$rho),
    paste0(nodes$display_name, "\n", nodes$name),
    paste0(nodes$display_name, "\n", nodes$name, "\nρ = ", sprintf("%.3f", nodes$rho), "\n(block ", nodes$block, ")")
  )
)

nodes$depth <- sapply(strsplit(nodes$name, "\\."), length)

# GRAPH AND LAYOUT
g <- graph_from_data_frame(edges_full, directed = TRUE, vertices = nodes)

layout_df <- create_layout(g, layout = "tree")

layout_df$true_depth <- nodes$depth[match(layout_df$name, nodes$name)]
layout_df$y <- layout_df$true_depth
layout_df <- as.data.frame(layout_df)


# LABEL FUNCTION FOR EXTRA NODES
make_label <- function(v) {
  
  if(v %in% c("XEC", "XDY")) {
    return(paste0(v, "\nrecombinant", "\nρ = ", sprintf("%.3f", rho_map[v]), "\n(block ", block_map[v], ")"))
  }
  
  if(is.na(rho_map[v])) {
    return(v)
  }
  
  paste0(v, "\n", unalias[v], "\nρ = ", sprintf("%.3f", rho_map[v]), "\n(block ", block_map[v], ")")
}

# ADD RECOMBINANTS, GHOST PARENTS, AND DUMMY PARENTS
# Recombinants are placed to the RIGHT of parents after coord_flip()

make_node <- function(template, name, label, x_shift, y_shift, recombinant = FALSE, ghost = FALSE, dummy = FALSE) {
  
  out <- template
  out$name <- name
  out$display_name <- name
  out$label <- label
  out$dummy <- dummy
  out$recombinant <- recombinant
  out$ghost <- ghost
  out$rho <- rho_map[name]
  out$block <- block_map[name]
  out$x <- template$x + x_shift
  out$y <- template$y + y_shift
  
  out
}

kp33_row <- layout_df %>% filter(display_name == "KP.3.3") %>% slice(1)
kp31_row <- layout_df %>% filter(display_name == "KP.3.1") %>% slice(1)

# XEC parents and recombinant

ks_row <- make_node(kp33_row, "KS.1.1", "KS.1.1", 1.5, 1.0, ghost = TRUE)
xec_row <- make_node(kp33_row, "XEC", make_label("XEC"), -0.5, 2.8, recombinant = TRUE)

d1_xec <- make_node(kp33_row, "Dummy1_XEC", "", -1.2, 1.6, dummy = TRUE)
d2_xec <- make_node(kp33_row, "Dummy2_XEC", "", -1.8, 1.2, dummy = TRUE)

# XDY parents and recombinant

kp32_row <- make_node(kp31_row, "KP.3.2", "KP.3.2", 4.0, 0.0, ghost = TRUE)
lb_row <- make_node(kp31_row, "LB.1.2.1", "LB.1.2.1", 3.0, 0.0, ghost = TRUE)

xdy_row <- make_node(kp31_row, "XDY", make_label("XDY"), 3.0, 1.5, recombinant = TRUE)

d1_xdy <- make_node(kp31_row, "Dummy1_XDY", "", 0.0, 1.7, dummy = TRUE)
d2_xdy <- make_node(kp31_row, "Dummy2_XDY", "", 0.0, 1.2, dummy = TRUE)

layout_df <- bind_rows(layout_df, ks_row, xec_row, d1_xec, d2_xec, kp32_row, lb_row, xdy_row, d1_xdy, d2_xdy)


# EDGE COORDINATES

main_edges <- edges_full %>%
  left_join(layout_df %>% select(name, x_from = x, y_from = y), by = c("from" = "name")) %>%
  left_join(layout_df %>% select(name, x_to = x, y_to = y), by = c("to" = "name"))

dummy_overlay <- dummy_edges %>%
  left_join(layout_df %>% select(display_name, x_from = x, y_from = y), by = c("from" = "display_name")) %>%
  left_join(layout_df %>% select(display_name, x_to = x, y_to = y), by = c("to" = "display_name"))

recomb_overlay <- recomb_edges %>%
  left_join(layout_df %>% select(display_name, x_from = x, y_from = y), by = c("from" = "display_name")) %>%
  left_join(layout_df %>% select(display_name, x_to = x, y_to = y), by = c("to" = "display_name"))


# SCALE BAR COORDINATES

scale_x <- min(layout_df$x, na.rm = TRUE) - 1.7
scale_y_start <- max(layout_df$y, na.rm = TRUE) - 2.0
scale_y_end <- scale_y_start + 1.0
scale_y_mid <- (scale_y_start + scale_y_end) / 2

# PLOT

set.seed(123)

p <- ggplot() +
  
  geom_segment(
    data = main_edges,
    aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
    linewidth = 0.8,
    color = "black"
  ) +
  
  geom_segment(
    data = dummy_overlay,
    aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
    linewidth = 0.8,
    color = "black",
    alpha = 0
  ) +
  
  geom_curve(
    data = recomb_overlay,
    aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
    curvature = 0.25,
    linetype = "dashed",
    linewidth = 1.1,
    color = "black",
    alpha = 1,
    arrow = arrow(length = unit(0.18, "cm"), type = "closed")
  ) +
  
  geom_point(
    data = layout_df %>% filter(!recombinant, !ghost),
    aes(x = x, y = y, fill = rho, alpha = dummy),
    shape = 21,
    size = 10,
    color = "black",
    stroke = 1
  ) +
  
  geom_point(
    data = layout_df %>% filter(ghost),
    aes(x = x, y = y),
    shape = 22,
    size = 11,
    fill = "grey30",
    color = "black",
    stroke = 1
  ) +
  
  geom_point(
    data = layout_df %>% filter(recombinant),
    aes(x = x, y = y, fill = rho),
    shape = 24,
    size = 9,
    color = "black",
    stroke = 1.2
  ) +
  
  geom_label_repel(
    data = layout_df,
    aes(x = x, y = y, label = label),
    size = 6.7,
    fontface = "bold",
    lineheight = 0.9,
    fill = "white",
    alpha = 1,
    label.size = 0.5,
    label.padding = unit(0.20, "lines"),
    label.r = unit(0.12, "lines"),
    box.padding = 1.0,
    point.padding = 0.7,
    segment.color = "grey50",
    segment.alpha = 0.6,
    max.overlaps = Inf
  ) +
  
# PANGO DISTANCE SCALE BAR


annotate(
  "segment",
  x = scale_x,
  xend = scale_x,
  y = scale_y_start,
  yend = scale_y_end,
  linewidth = 1.3,
  color = "white"
) +
  
  annotate(
    "segment",
    x = scale_x - 0.15,
    xend = scale_x + 0.15,
    y = scale_y_start,
    yend = scale_y_start,
    linewidth = 1.3,
    color = "white"
  ) +
  
  annotate(
    "segment",
    x = scale_x - 0.15,
    xend = scale_x + 0.15,
    y = scale_y_end,
    yend = scale_y_end,
    linewidth = 1.3,
    color = "white"
  ) +
  
  annotate(
    "text",
    x = scale_x - 0.45,
    y = scale_y_mid,
    label = "Genetic distance = 1",
    size = 9,
    fontface = "bold",
    color = "white"
  ) +
  
  scale_fill_gradient2(
    low = "darkblue",
    mid = "white",
    high = "red",
    midpoint = 1,
    na.value = "grey40",
    name = expression(rho)
  ) +
  
  scale_alpha_manual(values = c("TRUE" = 0, "FALSE" = 1), guide = "none") +
  
  coord_flip() +
  theme_void() +
  
  guides(
    fill = guide_colorbar(
      title.position = "top",
      barheight = unit(6, "cm"),
      barwidth = unit(0.8, "cm")
    )
  ) +
  
  theme(
    legend.position = c(0.99, 0.9),
    legend.title = element_text(size = 25, face = "bold", color = "black"),
    legend.text = element_text(size = 19, color = "black"),
    plot.margin = margin(20, 80, 20, 80)
  )

ggsave("lineage_tree.tiff",
  plot = p, device = "tiff", width = 25, height = 19, dpi = 600,bg = "white")
