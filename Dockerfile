# # Imagen base de shiny
FROM rocker/shiny-verse:latest

# # Instalando paquetes de ubuntu
RUN apt install software-properties-common -y; \
    add-apt-repository ppa:deadsnakes/ppa -y; \
    apt-get update && apt-get install -y \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    curl \
    g++ \
    build-essential \
    python3.11 \
    python3 \
    python3-dev \
    python3-distutils \
    python3-venv \
    python3-lib2to3 \
    python3-pip

# Preparando python

RUN python3 -m pip install --upgrade pip
RUN pip3 install virtualenv

# # Estableciendo el ambientes virtualese instalando dependencias
# Al parecer renv no funciona bien con shiny server, es necesario instalar las versiones manualmente
# Si existen problemas de compatibilidad con algun paquete se pueden consultar las versiones en renv.lock
RUN R -e "install.packages('renv', dependencies = TRUE, repos='https://cloud.r-project.org/')" \
    R -e "install.packages(c('arules', 'base64enc', 'bslib', 'cachem', 'caret', 'cli', 'clock'))" \
	R -e "install.packages(c('colorspace', 'commonmark', 'cpp11', 'crayon', 'data.table', 'digest'))" \
	R -e "install.packages(c('dplyr', 'e1071', 'ellipsis', 'fansi', 'farver', 'fastmap'))" \
	R -e "install.packages(c('fontawesome', 'foreach', 'fs', 'future', 'future.apply', 'gbm'))" \
	R -e "install.packages(c('generics', 'ggplot2', 'globals', 'glue', 'gower', 'gtable', 'hardhat'))" \
	R -e "install.packages(c('here', 'htmltools', 'httpuv', 'inTrees', 'ipred', 'isoband'))"  \
	R -e "install.packages(c('iterators', 'jquerylib', 'jsonlite', 'labeling', 'later', 'lava'))"  \
	R -e "install.packages(c('lifecycle', 'listenv', 'lubridate', 'magrittr', 'memoise', 'mime'))"  \
	R -e "install.packages(c('ModelMetrics', 'munsell', 'numDeriv', 'parallelly', 'pillar'))"  \
	R -e "install.packages(c('pkgconfig', 'plyr', 'png', 'pROC', 'prodlim', 'progressr', 'promises'))"  \
	R -e "install.packages(c('proxy', 'purrr', 'R6', 'randomForest', 'rappdirs', 'RColorBrewer'))"  \
	R -e "install.packages(c('Rcpp', 'RcppTOML', 'recipes', 'reshape2', 'reticulate', 'rlang'))"  \
	R -e "install.packages(c('rprojroot', 'RRF', 'RSNNS', 'sass', 'scales', 'shiny', 'shinyFeedback'))"  \
	R -e "install.packages(c('shinythemes', 'snn', 'sourcetools', 'SQUAREM', 'stringi', 'stringr'))"  \
	R -e "install.packages(c('tibble', 'tidyr', 'tidyselect', 'timechange', 'timeDate', 'tzdb'))"  \ 
	R -e "install.packages(c('utf8', 'vctrs', 'viridisLite', 'withr', 'wsrf', 'xgboost', 'xtable'))"

COPY ["requirements.txt", "/srv/shiny-server/webapp-shiny-mlbp/"]


WORKDIR /srv/shiny-server/webapp-shiny-mlbp

RUN pip install -r requirements.txt

# Abriendo puerto
EXPOSE 3838

# Activando shiny
CMD ["/usr/bin/shiny-server"]
