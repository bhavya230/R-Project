# 01_clean_data.R --------------------------------------------------------
# Phase 1a: import + clean  ->  data/clean_indonesia.rds
source("R/Setup.R")

raw_path <- "data/Dataset of AI Adoption Usage among Students in Indonesia Higher Education (1).xlsx"
raw <- read_excel(raw_path)

# 1. The sheet has summary tables pasted BELOW the respondents.
#    Respondents are rows 1:535 only.
dat <- raw %>%
  slice(1:535) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))

stopifnot(nrow(dat) == 535, max(dat$Respondent, na.rm = TRUE) == 535)

dat <- dat %>%
  rename(id = Respondent, gender = Gender, university = University,
         edu_level = `Level of Education`, province = Province,
         field = `Fields of Study`, ai_tool = `Type of AI`)

stopifnot(all(items %in% names(dat)), length(items) == 39)

# 2. The questionnaire is 4-point. Anything outside 1:4 is a data-entry error
#    (the raw file has three 5s: ATU3 x1, ATU5 x2) -> set to NA.
oor <- dat %>%
  summarise(across(all_of(items), ~ sum(!is.na(.x) & !(.x %in% 1:4)))) %>%
  pivot_longer(everything(), names_to = "item", values_to = "n_out_of_range") %>%
  filter(n_out_of_range > 0)
print(oor)

dat <- dat %>%
  mutate(across(all_of(items), ~ if_else(.x %in% 1:4, .x, NA_real_)))

# 3. Recodes ---------------------------------------------------------------
# Education: paper reports D3 = 91 and S1 = 444. The file has 435 coded 2 and
# 9 coded 3, and 435 + 9 = 444, so code 3 is treated as S1 (a miscode).
# CONFIRM against the Mendeley codebook.
dat <- dat %>%
  mutate(
    edu_level = if_else(edu_level == 3, 2, edu_level),
    female    = as.numeric(gender == 2),      # 1 = male, 2 = female (369 female)
    private   = as.numeric(university == 2),  # 1 = public, 2 = private (61)
    s1        = as.numeric(edu_level == 2),   # 0 = D3, 1 = S1
    chatgpt   = as.numeric(ai_tool == 3)
  )

stopifnot(sum(dat$female) == 369, sum(dat$private) == 61, sum(dat$s1) == 444)

# Field of study: labels INFERRED by matching counts to Table 1 of the paper
# (Social sciences = codes 0 + 2 = 121 + 76 = 197). Verify with the codebook.
dat <- dat %>%
  mutate(
    field_lab = case_when(
      field == 1 ~ "Natural sciences",
      field == 3 ~ "Engineering",
      field == 4 ~ "Mathematics",
      field == 5 ~ "Language",
      field == 6 ~ "Economics",
      field == 7 ~ "Education",
      field == 8 ~ "Art",
      field %in% c(0, 2) ~ "Social sciences"
    ),
    # Collapsed groups for Phase 3 (the raw fields include n = 11, 21, 30)
    field_grp = case_when(
      field_lab %in% c("Natural sciences", "Engineering", "Mathematics") ~ "STEM",
      field_lab %in% c("Social sciences", "Education")                   ~ "Social sci & education",
      TRUE                                                               ~ "Business, language & arts"
    )
  )

print(count(dat, field_lab))
print(count(dat, field_grp))

# AI tool labels, also inferred from counts in Table 1
dat <- dat %>%
  mutate(ai_tool_lab = recode(as.character(ai_tool),
                              "1" = "Dall-e", "2" = "Canva AI", "3" = "ChatGPT", "4" = "DeepL",
                              "5" = "Bing", "6" = "Tome", "7" = "Gamma", "9" = "Eduten",
                              "10" = "Shutterstock", "11" = "Perplexity", "12" = "Other"))

# 4. Missingness -----------------------------------------------------------
cat("\nMissing item cells:", sum(is.na(dat[items])),
    "of", nrow(dat) * length(items), "\n")
cat("Respondents with >=1 missing item:",
    sum(!complete.cases(dat[items])), "\n")
ggsave("output/missingness.png", gg_miss_var(dat[items]), width = 7, height = 8)

# 5. Descriptives ----------------------------------------------------------
desc <- describe(dat[items]) %>% as_tibble(rownames = "item")
write_csv(desc, "output/item_descriptives.csv")
print(desc %>% select(item, mean, sd, skew, kurtosis), n = 39)

# NOTE: if the questionnaire shows reverse-worded items (AUP1 looks suspicious:
# mean 1.99 vs 2.6-3.2 for its siblings), reverse them HERE and rerun:
# dat <- dat %>% mutate(AUP1 = 5 - AUP1)

saveRDS(dat, "data/clean_indonesia.rds")
