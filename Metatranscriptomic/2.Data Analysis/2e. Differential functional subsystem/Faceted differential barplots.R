##FIGURE 6 — Faceted differential barplots

# Load libraries
library(tidyverse)
library(stringr)

# -----------------------------
# 1. Read DESeq2 subsystem files
# -----------------------------

start  <- read.delim("start_level1.txt")
middle <- read.delim("middle_level1.txt")
end    <- read.delim("end_level1.txt")

# Add Stage column
start$Stage  <- "Start"
middle$Stage <- "Middle"
end$Stage    <- "End"

# Combine all stages
res_all <- bind_rows(start, middle, end)

# -----------------------------
# 2. Clean and annotate data
# -----------------------------

res_all <- res_all %>%
  filter(Row.names != "NO HIERARCHY") %>%
  mutate(
    Direction = ifelse(log2FoldChange > 0,
                       "Higher in Inoculated",
                       "Higher in Spontaneous"),
    Significant = ifelse(padj < 0.05,
                         "Significant",
                         "Not Significant")
  )

# Order Stage factor
res_all$Stage <- factor(res_all$Stage,
                        levels = c("Start", "Middle", "End"))

# -----------------------------
# 3. Keep subsystems with meaningful variation
#    (max |log2FC| > 0.5 in at least one stage)
# -----------------------------

res_filtered <- res_all %>%
  group_by(Row.names) %>%
  filter(max(abs(log2FoldChange)) > 0.5) %>%
  ungroup()

# -----------------------------
# 4. Order pathways based on Start stage
#    (biologically structured ordering)
# -----------------------------

order_df <- res_filtered %>%
  filter(Stage == "Start") %>%
  arrange(log2FoldChange)

ordered_pathways <- order_df$Row.names

res_filtered$Row.names <- factor(res_filtered$Row.names,
                                 levels = ordered_pathways)


# -----------------------------
# 5. Plot: Faceted differential barplot
# -----------------------------

ggplot(res_filtered,
       aes(x = Row.names,
           y = log2FoldChange,
           fill = Direction,
           alpha = Significant)) +
  geom_col(color = "black") +
  coord_flip() +
  facet_wrap(~ Stage) +
  scale_fill_manual(values = c(
    "Higher in Inoculated" = "#b2182b",
    "Higher in Spontaneous" = "#2166ac"
  )) +
  scale_alpha_manual(values = c(
    "Significant" = 1,
    "Not Significant" = 0.6
  )) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "black") +
  theme_bw() +
  labs(x = "",
       y = "log2 Fold Change (Inoculated vs Spontaneous)",
       fill = "",
       alpha = "") +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "top",
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 9)
  )

###Fixing some issues

library(tidyverse)

# Read files
start  <- read.delim("start_level1.txt")
middle <- read.delim("middle_level1.txt")
end    <- read.delim("end_level1.txt")

# Add stage column
start$Stage  <- "Start"
middle$Stage <- "Middle"
end$Stage    <- "End"

# Combine
res_all <- bind_rows(start, middle, end)

res_all <- res_all %>%
  filter(Row.names != "NO HIERARCHY") %>%
  mutate(
    Direction = ifelse(log2FoldChange > 0,
                       "Higher in Inoculated",
                       "Higher in Spontaneous"),
    Significant = ifelse(padj < 0.05, "Significant", "Not Significant")
  )


res_all$Stage <- factor(res_all$Stage,
                        levels = c("Start", "Middle", "End"))


##Keep subsystems that vary meaningfully
res_filtered <- res_all %>%
  group_by(Row.names) %>%
  filter(max(abs(log2FoldChange)) > 0.5)


ggplot(res_filtered,
       aes(x = reorder(Row.names, log2FoldChange),
           y = log2FoldChange,
           fill = Direction,
           alpha = Significant)) +
  geom_col(color = "black") +
  coord_flip() +
  facet_wrap(~ Stage, scales = "free_y") +
  scale_fill_manual(values = c(
    "Higher in Inoculated" = "#b2182b",
    "Higher in Spontaneous" = "#2166ac"
  )) +
  scale_alpha_manual(values = c(
    "Significant" = 1,
    "Not Significant" = 0.6
  )) +
  theme_bw() +
  labs(x = "",
       y = "log2 Fold Change (Inoculated vs Spontaneous)") +
  theme(
    legend.title = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )


##STEP 1 — Create global ordering

# Create ordering variable based on maximum absolute change across stages
order_df <- res_filtered %>%
  group_by(Row.names) %>%
  summarise(max_absFC = max(abs(log2FoldChange))) %>%
  arrange(max_absFC)

# Extract ordered names
ordered_pathways <- order_df$Row.names

# Apply global factor ordering
res_filtered$Row.names <- factor(res_filtered$Row.names,
                                 levels = ordered_pathways)

ggplot(res_filtered,
       aes(x = Row.names,
           y = log2FoldChange,
           fill = Direction,
           alpha = Significant)) +
  geom_col(color = "black") +
  coord_flip() +
  facet_wrap(~ Stage) +
  scale_fill_manual(values = c(
    "Higher in Inoculated" = "#b2182b",
    "Higher in Spontaneous" = "#2166ac"
  )) +
  scale_alpha_manual(values = c(
    "Significant" = 1,
    "Not Significant" = 0.6
  )) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "",
       y = "log2 Fold Change (Inoculated vs Spontaneous)") +
  theme(
    legend.title = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )


##PART 2 — Increase readability of subsystem names

library(stringr)

res_filtered$Row.names <- str_wrap(res_filtered$Row.names, width = 30)

res_filtered$Row.names <- recode(res_filtered$Row.names,
                                 "Cofactors, Vitamins, Prosthetic Groups, Pigments" =
                                   "Cofactors & Vitamins",
                                 "Phages, Prophages, Transposable elements, Plasmids" =
                                   "Phages & Mobile Elements"
)



order_df <- res_filtered %>%
  group_by(Row.names) %>%
  summarise(meanFC = mean(log2FoldChange)) %>%
  arrange(meanFC)

ordered_pathways <- order_df$Row.names

res_filtered$Row.names <- factor(res_filtered$Row.names,
                                 levels = ordered_pathways)







##More simple 
library(tidyverse)

# Read files
start  <- read.delim("start_level1.txt")
middle <- read.delim("middle_level1.txt")
end    <- read.delim("end_level1.txt")

# Add stage column
start$Stage  <- "Start"
middle$Stage <- "Middle"
end$Stage    <- "End"

# Combine
res_all <- bind_rows(start, middle, end)


res_sig <- res_all %>%
  filter(pvalue < 0.5)

res_sig <- res_sig %>%
  filter(Row.names != "NO HIERARCHY")


ggplot(res_sig,
       aes(x = reorder(Row.names, log2FoldChange),
           y = log2FoldChange,
           fill = log2FoldChange > 0)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ Stage, scales = "free_y") +
  scale_fill_manual(values = c("TRUE" = "#b2182b",
                               "FALSE" = "#2166ac"),
                    labels = c("FALSE" = "Higher in Spontaneous",
                               "TRUE"  = "Higher in Inoculated")) +
  theme_bw() +
  labs(x = "",
       y = "log2 Fold Change (Inoculated vs Spontaneous)",
       fill = "")



table(res_all$Stage)

table(res_sig$Stage)

res_all <- res_all %>%
  filter(Row.names != "NO HIERARCHY") %>%
  mutate(Significant = ifelse(padj < 0.05, "Yes", "No"))

ggplot(res_all,
       aes(x = reorder(Row.names, log2FoldChange),
           y = log2FoldChange,
           fill = Significant)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ Stage, scales = "free_y") +
  scale_fill_manual(values = c("Yes" = "#b2182b",
                               "No" = "grey70")) +
  theme_bw() +
  labs(x = "",
       y = "log2 Fold Change (Inoculated vs Spontaneous)")


table(res_sig$Stage)


ggplot(res_all,
       aes(x = reorder(Row.names, log2FoldChange),
           y = log2FoldChange,
           fill = log2FoldChange > 0)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ Stage, scales = "free_y") +
  scale_fill_manual(values = c("TRUE" = "#b2182b",
                               "FALSE" = "#2166ac"),
                    labels = c("FALSE" = "Higher in Spontaneous",
                               "TRUE"  = "Higher in Inoculated")) +
  theme_bw() +
  labs(x = "",
       y = "log2 Fold Change (Inoculated vs Spontaneous)",
       fill = "")


##FIGURE 2 — Global heatmap of functional shifts
library(pheatmap)

# Keep all subsystems (even non-significant)
res_heat <- res_all %>%
  filter(Row.names != "NO HIERARCHY") %>%
  select(Row.names, Stage, log2FoldChange)

# Convert to wide format
heat_matrix <- res_heat %>%
  pivot_wider(names_from = Stage,
              values_from = log2FoldChange)

# Set rownames
heat_matrix <- as.data.frame(heat_matrix)
rownames(heat_matrix) <- heat_matrix$Row.names
heat_matrix$Row.names <- NULL

# Convert to matrix
heat_matrix <- as.matrix(heat_matrix)

pheatmap(heat_matrix,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main = "Functional pathway shifts across fermentation stages")


