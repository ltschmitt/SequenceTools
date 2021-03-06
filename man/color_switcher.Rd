% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/color_switcher.R
\name{color_switcher}
\alias{color_switcher}
\title{Switch the order of the colors}
\usage{
color_switcher(n, color_function = gg_color_hue, ordering = "systematic")
}
\arguments{
\item{n}{integer. number of colors to produce}

\item{color_function}{function to use for generating colors}

\item{ordering}{'systematic' or 'random'. the method it switches the colors}
}
\value{
character vector with color values
}
\description{
This function takes a color function and reorders the colors that it outputs.
}
\examples{
color_switcher(6, ordering = 'random')
}
\keyword{color,}
\keyword{plotting}
