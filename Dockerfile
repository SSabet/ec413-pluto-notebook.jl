FROM jupyter/base-notebook:latest

USER root

RUN pip install --no-cache --upgrade pip

RUN apt-get update && \
    apt-get install -y --no-install-recommends wget && \
    apt-get install -y --no-install-recommends build-essential && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

#USER ${NB_USER}
ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

# RUN adduser --disabled-password \
#     --gecos "Default user" \
#     --uid ${NB_UID} \
#     ${NB_USER}
USER ${NB_USER}

COPY --chown=${NB_USER}:users ./plutoserver ./plutoserver
COPY --chown=${NB_USER}:users ./environment.yml ./environment.yml
COPY --chown=${NB_USER}:users ./setup.py ./setup.py
COPY --chown=${NB_USER}:users ./runpluto.sh ./runpluto.sh
COPY --chown=${NB_USER}:users ./notebooks ./notebooks
COPY --chown=${NB_USER}:users ./images ./images
COPY --chown=${NB_USER}:users ./Project.toml ./Project.toml
COPY --chown=${NB_USER}:users ./Manifest.toml ./Manifest.toml
COPY --chown=${NB_USER}:users ./warmup.jl ./warmup.jl
COPY --chown=${NB_USER}:users ./create_sysimage.jl ./create_sysimage.jl

USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

RUN jupyter labextension install @jupyterhub/jupyter-server-proxy && \
    jupyter lab build && \
    jupyter lab clean && \
    pip install . --no-cache-dir && \
    #pip install matplotlib && \
    rm -rf ~/.cache

USER root

RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.9/julia-1.9.3-linux-x86_64.tar.gz && \
    tar -xvzf julia-1.9.3-linux-x86_64.tar.gz && \
    mv julia-1.9.3 /opt/ && \
    ln -s /opt/julia-1.9.3/bin/julia /usr/local/bin/julia && \
    rm julia-1.9.3-linux-x86_64.tar.gz

RUN apt-get update && \
    apt-get install -y texlive-latex-extra dvipng && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

USER ${NB_USER}

ENV USER_HOME_DIR /home/${NB_USER}
ENV JULIA_PROJECT ${USER_HOME_DIR}
ENV JULIA_DEPOT_PATH ${USER_HOME_DIR}/.julia
WORKDIR ${USER_HOME_DIR}

RUN julia --project=${USER_HOME_DIR} -e "import Pkg; Pkg.Registry.update(); Pkg.instantiate(); Pkg.precompile()"

RUN julia --project=${USER_HOME_DIR} create_sysimage.jl
RUN julia -J${USER_HOME_DIR}/sysimage.so --project=${USER_HOME_DIR} -e "import Pkg; Pkg.precompile()"
RUN julia --project=${USER_HOME_DIR} -e "import Pkg; Pkg.precompile()"

