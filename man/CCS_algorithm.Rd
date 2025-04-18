% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CCS_algorithm.R
\name{CCS_algorithm}
\alias{CCS_algorithm}
\title{Cumulative Continuity Score (CCS) algorithm}
\usage{
CCS_algorithm(x, value_col, threshold, full_result)
}
\arguments{
\item{x}{a data frame or tibble.}

\item{value_col}{character; the column in which the change points will be analysed. This column should contain binary data, e.g. 0 and 1.}

\item{threshold}{numeric; the CCS value above which the segments will be included for locating the change points. Default = 50.}

\item{full_result}{logical; it indicates whether the result having CCS at all rows will be produced, or only the change points otherwise.}
}
\description{
Cumulative Continuity Score (CCS) algorithm for locating change points.
}
