\name{plot_replicates}
\alias{plot_replicates}
\title{Show the ensemble fishstock and harvest dynamics}
\usage{
  plot_replicates(sims, S = NULL, allee = 0, e_star = NULL,
    fish_reps = TRUE, fish_sd = TRUE, fish_mean = TRUE,
    harvest = TRUE, harvest_reps = TRUE,
    harvest_mean = TRUE, harvest_sd = TRUE,
    crash_count = TRUE)
}
\arguments{
  \item{sims}{simulation results}

  \item{S}{Reed's S, from opt$S, where opt <-
  find_dp_optim(...)}

  \item{allee}{the allee threshold (defaults to 0)}

  \item{fish_reps}{logical - show replicate trajectories of
  fish stocks?}

  \item{fish_sd}{logical - show ribbon of SD of fishstock?}

  \item{fish_mean}{logical - show mean fishstock?}

  \item{harvest}{- include harvest level on the same plot?}

  \item{harvest_reps}{- show individual reps for harvest?}

  \item{harvest_mean}{- show mean harvest level over time?}

  \item{harvest_sd}{- show harvest SD}

  \item{crash_count}{- count the number of reps that have
  fallen below allee}
}
\value{
  a ggplot object (display with print(object))
}
\description{
  Show the ensemble fishstock and harvest dynamics
}
\details{
  if line values for any S, allee, e_star are not specified
  then these are omitted from the plot.
}

