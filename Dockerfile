FROM mavelli/rocker-bayesian

WORKDIR /app

RUN R -e "install.packages(c( \
  'boot', \
  'doParallel', \
  'dplyr', \
  'fdrtool', \
  'ggplot2', \
  'here', \
  'jagsUI', \
  'rjags', \
  'pracma', \
  'readxl', \
  'reshape2', \
  'scales', \
  'stringr', \
  'tictoc', \
  'tidyr', \
  'vctrs', \
  dependencies=TRUE \
))"