% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prodd_parameters.R
\name{transform_parameters}
\alias{transform_parameters}
\title{Update the experimental_design of a parameter object}
\usage{
transform_parameters(params, new_experimental_design)
}
\arguments{
\item{params}{an object of class `prodd_parameters` which is returned by the
\code{fit_parameters()} function.}

\item{new_experimental_design}{a vector that assignes each sample to one condition.
It has to have the same length as the old experimental_design. It can either be
a factor, a character or a numeric vector. Each unique element is one condition.
To remove samples code them as \code{NA}.}
}
\value{
the updated list of class prodd_parameters
}
\description{
Helper method to transform the parameters from one experimental_design
to another. It avoids expensive re-calculation of the parameters and
instead uses a averaging of the parameters to produce the new values.
}
\details{
The method only works for changing the condition assignment of samples
or removing samples. The method cannot add new samples and infer the
parameters for them.

The method is for example useful if you want to calculate the distance
with unbiased feature_parameters that don't already include the information
to which condition a sample belongs.
}
