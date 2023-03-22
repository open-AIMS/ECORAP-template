## ---- setup
assignInNamespace('.sep.label',  "^\\ *(#|--)+\\s*(@knitr|----+)(.*?)-*\\s*$", ns='knitr')

tidyverse_style_with_comments_removed <- function() {
  remove_comments <- function(pd) {
    is_comment <- pd$token == "COMMENT"
    pd <- pd[!is_comment,]
    pd
  }
  tidyverse_with_comments_removed <- styler::tidyverse_style()
  tidyverse_with_comments_removed$token$remove_comments <- remove_comments
  tidyverse_with_comments_removed
}

knitr::opts_chunk$set(echo = TRUE, message=FALSE, warning=FALSE,cache.lazy = FALSE, tidy='styler',
                      tidy.opts=list(transformers = tidyverse_style_with_comments_removed()))
## knitr::opts_chunk$set(echo = TRUE, warning=FALSE, message=FALSE)

# save the built-in output hook
hook_output <- knitr::knit_hooks$get("output")

# set a new output hook to truncate text output
knitr::knit_hooks$set(output = function(x, options) {
  if (!is.null(n <- options$out.lines)) {
    x <- xfun::split_lines(x)
    if (length(x) > n) {
      # truncate the output
      x <- c(head(x, n), "....\n")
    }
    x <- paste(x, collapse = "\n")
  }
  hook_output(x, options)
})
options(tinytex.engine = 'xelatex')

## ----end

## ---- loadPackages
library(knitr)
library(tidyverse)
library(simstudy)
library(ggraph)
library(igraph)
library(forcats)
library(patchwork)
library(ggbrace)

library(brms)
library(INLA)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(ggeffects)
library(posterior)
library(HDInterval)
library(tidybayes)
## ----end

## ---- preparePaths
DATA_PATH <<- "../data/"
OUTPUT_PATH <<- "../output/"
FIGS_PATH <<- paste0(OUTPUT_PATH, "figures")
TABS_PATH <<- paste0(OUTPUT_PATH, "tables")

if (!dir.exists(DATA_PATH)) dir.create(DATA_PATH)
if (!dir.exists(paste0(DATA_PATH,"primary"))) dir.create(paste0(DATA_PATH, "primary"))
if (!dir.exists(paste0(DATA_PATH,"processed"))) dir.create(paste0(DATA_PATH, "processed"))
if (!dir.exists(paste0(DATA_PATH,"modelled"))) dir.create(paste0(DATA_PATH, "modelled"))
if (!dir.exists(paste0(DATA_PATH,"summarised"))) dir.create(paste0(DATA_PATH, "summarised"))
if (!dir.exists(paste0(DATA_PATH,"workspace"))) dir.create(paste0(DATA_PATH, "workspace"))

if (!dir.exists(OUTPUT_PATH)) dir.create(OUTPUT_PATH)
if (!dir.exists(FIGS_PATH)) dir.create(FIGS_PATH)
if (!dir.exists(TABS_PATH)) dir.create(TABS_PATH)
## ----end
