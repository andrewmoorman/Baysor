FROM debian:bullseye

ARG BAYSOR_VER=v0.6.2
LABEL baysor=$BAYSOR_VER
ARG JULIA_VER=1.8.3
LABEL julia=$JULIA_VER

RUN apt-get update && apt-get install -y build-essential zip unzip wget git
RUN mkdir /depot
ENV JULIA_PATH=/usr/local/julia
ENV JULIA_DEPOT_PATH=/depot
ENV JULIA_PROJECT=/baysor

### Install Julia and add to PATH
# NOTE: Only successfully tested with Julia <= v1.8.3
RUN wget \
https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_VER%.*}/\
julia-${JULIA_VER}-linux-x86_64.tar.gz -P ${JULIA_PATH}
RUN cd ${JULIA_PATH} && \
tar -xzf ${JULIA_PATH}/julia-${JULIA_VER}-linux-x86_64.tar.gz
ENV PATH="$PATH:${JULIA_PATH}/julia-${JULIA_VER}/bin"

### Compile Baysor
# NOTE: Only successfully tested with Baysor v0.5.0, v0.6.2
ARG REPO=https://github.com/kharchenkolab/Baysor.git
RUN mkdir ${JULIA_PROJECT} && cd ${JULIA_PROJECT} && \
git clone --branch ${BAYSOR_VER} ${REPO}
RUN julia --project=${JULIA_PROJECT}/Baysor -e '\
using Pkg; \
Pkg.instantiate(); \
Pkg.precompile(); \
Pkg.build()'

RUN chmod 777 -R /depot/

ENV PATH="${JULIA_DEPOT_PATH}/bin:${PATH}"

ENV JULIA_DEPOT_PATH="/tmp/:${JULIA_DEPOT_PATH}"

ENTRYPOINT ["/bin/bash"]