#!/bin/sh

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")"

version=4.4.0
workdir=$(pwd)/rstudio-container-tmp

mkdir -p -m 700 ${workdir}/${version}/run ${workdir}/${version}/tmp ${workdir}/${version}/var/lib/rstudio-server
cat > ${workdir}/database.conf <<END
provider=sqlite
directory=/var/lib/rstudio-server
END

# Do not suspend idle sessions.
# Alternative to setting session-timeout-minutes=0 in /etc/rstudio/rsession.conf
# https://github.com/rstudio/rstudio/blob/v1.4.1106/src/cpp/server/ServerSessionManager.cpp#L126
export APPTAINERENV_RSTUDIO_SESSION_TIMEOUT=0
export APPTAINERENV_USER=$(id -un)

singularity exec --cleanenv -c -W ${workdir}/${version} \
                 --bind .:/home/$USER,/research:/research,/hpcf:/hpcf,${workdir}/${version}/run:/run,${workdir}/${version}/tmp:/tmp,${workdir}/database.conf:/etc/rstudio/database.conf,${workdir}/${version}/var/lib/rstudio-server:/var/lib/rstudio-server rstudio_4.4.0_seurat_4.4.0_latest.sif \
  /bin/bash
