htmlTable <- function(
  df,
  caption="Basic Infos",
  col.names = NA,
  row.names = NA
){
  require(kableExtra)
  kable_styling(knitr::kable(df,
                             format = "html",
                             col.names = col.names,
                             row.names = row.names,
                             align = "c",
                             caption = caption,
                             table.attr = "class=\"table table-bordered\""),
                bootstrap_options = c("striped", "hover", "condensed", "responsive"))
}
