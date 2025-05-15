
# Rashnu

**Rashnu** is the Zoroastrian deity of truth and justice—the one who
weighs souls on a golden scale.  
This R package, developed by **Zarathu**, draws inspiration from
Rashnu’s role as the divine judge to offer precision and fairness in
sample size determination.

> *“Where truth is weighed, science begins.”*

------------------------------------------------------------------------

In clinical trials and research design, every decision matters.  
**Rashnu** helps researchers define the *right* number of participants
for: - Non-inferiority studies - Superiority comparisons (Lakatos
method) - One-arm survival designs with transformation-based inference

This package brings clarity, rigor, and justice to your design process.

## Installation

You can install the development version of rashnu from
[GitHub](https://github.com/zarathucorp/rashnu) with:

``` r
# install.packages("pak")
pak::pak("zarathucorp/rashnu")
```

## Example

### Rstudio Addins

Use `rashnuBasic()` or Click **Interactive Sample Size Calculator**
Addin

``` r
rashnuBasic()
```

![](inst/www/figures/addin.gif)

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

![](README_files/figure-gfm/pressure-1.png)<!-- -->

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
