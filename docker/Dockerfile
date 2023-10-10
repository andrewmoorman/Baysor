FROM debian:bullseye

RUN apt-get update && apt-get install -y build-essential zip unzip wget

### Install Julia and add to PATH
# NOTE: Only successfully tested with Julia <= v1.8.3
ARG JULIA_VERSION=1.8.3
RUN wget \
https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_VERSION%.*}/\
julia-${JULIA_VERSION}-linux-x86_64.tar.gz -P /root/
RUN cd root; tar -xzf /root/julia-${JULIA_VERSION}-linux-x86_64.tar.gz
ENV PATH="$PATH:/root/julia-${JULIA_VERSION}/bin"

# ### Download and Compile Baysor from zip
# # Download GitHub zip file
# ARG Baysor='https://github.com/kharchenkolab/Baysor/archive/refs/heads/master.zip'
# RUN wget ${Baysor}; unzip master.zip -d /root/; cd /root/Baysor-master
# # Compile Baysor
# RUN cd /root/Baysor-master; julia -e '\
# using Pkg; \
# Pkg.develop(path="/root/Baysor-master/"); \
# Pkg.build()'

### Compile Baysor
# NOTE: Only successfully tested with Baysor v0.5.0, v0.6.2
ARG BAYSOR_VERSION=v0.6.2
RUN julia -e "\
using Pkg; \
Pkg.add(\
url=\"https://github.com/kharchenkolab/Baysor.git\", \
rev=\"${BAYSOR_VERSION}\"\
); \
import Baysor; \
Pkg.activate(dirname(dirname(pathof(Baysor)))); \
Pkg.instantiate(); \
Pkg.build()"

# Add Baysor executable
RUN echo "export PATH=/root/.julia/bin/:$PATH" >> ~/.bashrc
RUN echo "alias julia='/usr/local/julia/bin/julia \
--sysimage=/root/.julia/scratchspaces/*/sysimg/libbaysor.so'" \
>> ~/.bashrc
RUN ln -s /root/.julia/bin/baysor /usr/local/bin/baysor

ENTRYPOINT ["/bin/bash"]
WORKDIR /root/