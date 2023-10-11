#' Generate G matrix from text input
#' 
#' Convenience function for parsing a multi-line text string into the G matrix, where each line 
#' specifies one causal relationship between traits. Each line must be formatted with four 
#' space-separated tokens as: 
#' 
#' <trait 1> <causal orientation> <trait 2> <effect>
#' 
#' The function tries to parse the effect from numeric expressions or global variables.
#' 
#' @param text A multi line text string that specifies the causal relationships between traits
#'
#' @export
#' @return A square matrix with dimension number of traits
#' 
#' @examples
#' 
#' # Simple specification
#' text <- "X <- Y 0.25
#'          A <- X 0.1
#'          A <- Y sqrt(0.34)"
#' G_from_text(text)
#' 
#' # Another example
#' a <- 0.25
#' b <- -0.34
#' text <- "X -> Y a
#'          A -> X b
#'          Y <- A sqrt(0.34)"
#' G_from_text(text)
G_from_text <- function(text) {
  f5f731c87529 <- strsplit(text, "\n")[[1]]
  f5f75f83f047 <- lapply(f5f731c87529, \(x) strsplit(trimws(x), " +")[[1]])
  if(!all(sapply(f5f75f83f047, length) == 4)) {
    stop("Expect there to be four tokens, space separated, for each line e.g. 'A <- B 0.25")
  }
  lapply(f5f75f83f047, \(x) {
    if(x[2] == "<-") {
      dplyr::tibble(eff = eval(parse(text=x[4])), i = x[3], j = x[1])
    } else if(x[2] == "->") {
      dplyr::tibble(eff = eval(parse(text=x[4])), i = x[1], j = x[3])
    } else {
      stop("Nodes must be linked by either '<-' or '->'. See Vignette for examples.")
    }
  }) %>% 
    dplyr::bind_rows() %>%
    G_from_df() %>%
    return()
}

#' Generate G matrix from data frame
#' 
#' Specify G matrix as a data frame in long format, one row represents one causal 
#' relationship of trait `i` on trait `j` with effect `eff`.
#' 
#' @param df A data frame with required columns i (the causal trait name), j (the 
#' response trait name) and eff (the effect of i on j)
#'
#' @export
#' @return A square matrix with dimension number of traits
#' 
#' @examples
#' 
#' # Simple example
#' df <- dplyr::tribble(~i, ~j, ~eff,
#'                      "Y", "X", 0.25,
#'                      "X", "A", 0.24,
#'                      "Y", "A", 0.34)
#' 
#' G_from_df(df)
#' 
G_from_df <- function(df) {
  stopifnot(inherits(df, "data.frame"))
  stopifnot(nrow(df) >0)
  stopifnot(all(c("i", "j", "eff") %in% names(df)))
  nodes <- unique(c(df$i, df$j))
  message(nrow(df), " causal relationships specified amongst ", length(nodes), " traits")
  G <- matrix(0, length(nodes), length(nodes))
  rownames(G) <- nodes
  colnames(G) <- nodes
  for(x in 1:nrow(df)) {
    G[df$i[x], df$j[x]] <- df$eff[x]
  }
  return(G)
}
