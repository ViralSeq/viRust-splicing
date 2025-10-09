library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
output_file <- args[2]

data <- read_tsv(input_file)

unspliced <- c(
  "D1_D1-unspliced",
  "D1_gag-AUG",
  "noD1_gag-AUG_unspliced"
)

# Function to get most common value (mode)
get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

data_with_family <- data %>% filter(umi_family != "None")

data_with_family <- data_with_family %>%
  group_by(umi_family, splice_category) %>%
  mutate(
    alternative_d1_used_common = get_mode(alternative_d1_used)
  ) %>%
  ungroup() %>%
  mutate(
    spliced = case_when(
      str_detect(splice_category, "D1_unknown") ~ "unknown",
      str_detect(splice_category, "no-acceptor") ~ "unknown",
      splice_category %in% unspliced ~ "unspliced",
      TRUE ~ "spliced"
    ),
    D1_type = case_when(
      (spliced == "spliced" & str_detect(splice_category, "noD1_")) ~
        "Alternative D1",
      (spliced == "spliced" & str_detect(splice_category, "D1prime_")) ~
        "D1prime",
      spliced == "spliced" ~ "D1"
    ),
    splice_category_2 = ifelse(
      spliced == "spliced",
      splice_category %>%
        str_replace("D1prime", "D1") %>%
        str_replace("noD1", "D1"),
      spliced
    )
  )


data_with_family_unique <-
  data_with_family %>%
  select(
    c(3:6, alternative_d1_used_common, spliced, D1_type, splice_category_2)
  ) %>%
  unique()

write_csv(
  data_with_family_unique,
  output_file
)