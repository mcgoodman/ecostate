% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ode23.R
\name{ode23}
\alias{ode23}
\title{Non-stiff (and stiff) ODE solvers}
\usage{
ode23(f, a, b, y0, n, Pars, rtol = 0.001, atol = 1e-06)
}
\arguments{
\item{f}{function in the differential equation \eqn{y' = f(x, y)};
defined as a function \eqn{R \times R^m \rightarrow R^m}, where \eqn{m} is the number of equations.}

\item{a}{starting time for the interval to integrate}

\item{b}{ending time for the interval to integrate.}

\item{y0}{starting values at time \code{a}}

\item{n}{Not used}

\item{Pars}{named list of parameters passed to f}

\item{rtol}{relative tolerance.}

\item{atol}{absolute tolerance.}
}
\value{
List with components t for time points between a and b and y an
n-by-m matrix with solutions for variables in columns, i.e.
each row contains one time stamp.
}
\description{
Runge-Kutta (2, 3)-method with variable step size, resp
}
\details{
Copied from pracma under GPL-3, with small modifications to work with RTMB.
This can be used to simulate dynamics, but not during estimation
}
