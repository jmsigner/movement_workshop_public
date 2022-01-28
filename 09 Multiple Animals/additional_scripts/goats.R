library(glmmTMB)
library(broom.mixed)
data(goats, package = "ResourceSelection")

# @knitr goat.models
# This is a naive approach (ignoring different animals)
m1 <- glmmTMB(STATUS ~ ELEVATION + SLOPE, 
              data = goats, family = binomial())

# This is the random intercept model
m2 <- glmmTMB(STATUS ~ ELEVATION + SLOPE + (1 | ID), 
              data = goats, family = binomial())

# This is a random slope and intercept model
m3 <- glmmTMB(STATUS ~ ELEVATION + SLOPE + 
                (ELEVATION + SLOPE | ID),
              data = goats, family = binomial())

# @knitr goaway
bind_rows(
  tidy(m1, conf.int = TRUE) %>% mutate(what = "glm"),
  tidy(m2, conf.int = TRUE) %>% mutate(what = "glmm (intercept)"),
  tidy(m3, conf.int = TRUE) %>% mutate(what = "glmm (intercept & slope)")
) %>% filter(effect == "fixed") %>% 
  mutate(term1 = str_remove_all(term, "[\\(\\)]"), 
         term1 = str_to_title(term1) %>% factor() %>% fct_inorder(), 
         what = factor(what) %>% fct_inorder()) %>% 
  ggplot(aes(what, estimate)) + 
           geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  facet_wrap(~ term1, scale = "free") + 
  geom_hline(yintercept = 0, col = "red", lty = 2) + 
  labs(y = "Estimate", x = "Model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank())

ggsave("img/rsf_re.png", width = 12, height = 8, units = "cm")
         