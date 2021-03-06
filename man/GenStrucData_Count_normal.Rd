% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count_normal_data.R
\name{GenStrucData_Count_normal}
\alias{GenStrucData_Count_normal}
\title{Generate count data}
\usage{
GenStrucData_Count_normal(u, v, n = 50, n_y = 1000, n_x = 1000,
  theta = 0.01, sigma_z = 1, sigma_x = NULL, sigma_y = NULL)
}
\arguments{
\item{u:}{caonical coefficients for x}

\item{v:}{canonical coefficients for y}

\item{n:}{sample size}

\item{sigma_z:}{sd of the latent variable z, scalar}

\item{sigma_x:}{covariance matrix of x, default to be identity}

\item{sigma_y:}{covariance matrix of y, default to be identity}

\item{n_y:}{parameter n in DM model}

\item{n_x:}{parameter n in DM model}

\item{theta:}{oversidpersion parameter in DM model}
}
\value{
a list containing data matrix X and Y
}
\description{
Generate count data, using normal data from latent variable model
}
