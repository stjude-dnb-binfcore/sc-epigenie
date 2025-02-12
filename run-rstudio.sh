#!/bin/sh

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")"

version=4.4.0
workdir=./rstudio-container-tmp

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

# Get unused socket per https://unix.stackexchange.com/a/132524
# Tiny race condition between the python & singularity commands
readonly PORT=$(python -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()')
# Get node IP address.
readonly ADD=$(nslookup `hostname` | grep -i address | awk -F" " '{print $2}' | awk -F# '{print $1}' | tail -n 1)

cat 1>&2 <<END
"Running RStudio at $ADD:$PORT"
END

# Singularity call to start RStudio Container
singularity exec --cleanenv -c -W ${workdir}/${version} \
                 --bind .:/home/$USER,/research:/research,/hpcf:/hpcf,${workdir}/${version}/run:/run,${workdir}/${version}/tmp:/tmp,${workdir}/database.conf:/etc/rstudio/database.conf,${workdir}/${version}/var/lib/rstudio-server:/var/lib/rstudio-server rstudio_4.4.0_seurat_4.4.0_latest.sif \
    rserver --www-port ${PORT} \
            --secure-cookie-key-file "/tmp/rstudio-server/secure-cookie-key" \
            --auth-stay-signed-in-days=30 \
            --auth-none=1 \
            --server-user $USER \
            --auth-timeout-minutes=0 
