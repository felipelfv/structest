# Confidence intervals for structest_fit parameters

Confidence intervals for structest_fit parameters

## Usage

``` r
# S3 method for class 'structest_fit'
confint(object, parm, level = object$conf_level, ...)
```

## Arguments

- object:

  an object of class `"structest_fit"`.

- parm:

  character vector of parameter names, or numeric indices. If missing,
  all parameters are included.

- level:

  confidence level (default uses the level from the fit).

- ...:

  further arguments (ignored).

## Value

A matrix with columns for lower and upper confidence bounds.
